### Main GP functions

#' gp_optimize
#'
#' to optimizes S2 using MLE. This function automatically runs when gp_train(..., optimize=TRUE).
#'
#' @param K kernel matrix (of covariates)
#' @param Y Y vector (outcome variable)
#'
#' @importFrom stats sd optimize
#' @importFrom Rcpp sourceCpp
#' @return \item{s2opt}{optimized s2 value}
#' \item{nlml}{the likelihood value at the optimal point}
#' @export
gp_optimize <- function(K, Y, optim.tol=0.1) {

  # Define the objective function
  nlml <- function(s2, K, Y) {
    return(-1 * log_marginal_likelihood_cpp(K, Y, s2))
  }

  # Find the optimal s2 (by MLE)
  # Since we always scale y, it is between 0 and 1; is 0.05 good for lower bound?
  opt <- optimize(nlml, interval=c(0.05, 1), K=K, Y=Y, maximum=FALSE, tol=optim.tol)
  results <- list(s2opt = opt$minimum,
                  nlml = opt$objective)
  return(results)
}


#' gp_train
#'
#' to train GP model with training data set
#'
#' @param X a set of covariate data frame or matrix
#' @param Y Y vector (outcome variable)
#' @param b bandwidth (default = NULL)
#' @param s2 noise or a fraction of Y not explained by X (default = 0.3)
#' @param optimize a logical value to indicate whether an automatic optimized value of S2 should be used. If FALSE, users must define s2. (default = FALSE)
#' @param scale a logical value to indicate whether covariates should be scaled. (dafault = TRUE)
#' @param kernel_type a character value indicating the kernel type (default = "gaussian")
#' @param period a numeric value for the period parameter, required for periodic kernels (default = NULL)
#' @param time_col a character or numeric value indicating which column is the time variable. If specified, this column will be moved to the first position for correct period scaling in periodic kernels. (default = NULL)
#' @param mixed_data a logical value to indicate whether the covariates contain a categorical/binary variable (default = FALSE)
#' @param cat_columns a character or a numerical vector indicating categorical variables. Must be character (not numeric) when time_col is specified. (default = NULL)
#' @param Xtest a data frame or a matrix of testing covariates. This is necessary when a non-overlapping categorical value exists between training and testing data sets. (default = NULL)
#'
#' @return \item{post_mean_scaled}{posterior distribution of Y in a scaled form}
#' \item{post_mean_orig}{posterior distribution of Y in an original scale}
#' \item{post_cov_scaled}{posterior covariance matrix in a scaled form}
#' \item{post_cov_orig}{posterior covariance matrix in an original scale}
#' \item{K}{a kernel matrix of X}
#' \item{prior_mean_scaled}{prior distribution of mean in a scaled form}
#' \item{X.orig}{the original matrix or data set of X}
#' \item{X.init}{the original matrix or data set of X with categorical variables in an expanded form}
#' \item{X.init.mean}{the initial mean values of X}
#' \item{X.init.sd}{the initial standard deviation values of X}
#' \item{Y.init.mean}{the initial mean value of Y}
#' \item{Y.init.sd}{the initial standard deviation value of Y}
#' \item{K}{the kernel matrix of X}
#' \item{Y}{scaled Y}
#' \item{X}{scaled X}
#' \item{b}{bandwidth}
#' \item{s2}{sigma squared}
#' \item{alpha}{alpha value in Rasmussen and Williams (2006) p.19}
#' \item{L}{L value in Rasmussen and Williams (2006) p.19}
#' \item{mixed_data}{a logical value indicating whether X contains a categorical/binary variable}
#' \item{cat_columns}{a character or a numerical vector indicating the location of categorical/binary variables in X}
#' \item{cat_num}{a numerical vector indicating the location of categorical/binary variables in an expanded version of X}
#' \item{time_col}{the time column specification used (or NULL if not specified)}
#' \item{Xcolnames}{column names of X}
#' @importFrom stats sd
#' @importFrom Rcpp sourceCpp
#'
#' @export
gp_train <- function(X, Y, b = NULL, s2 = 0.3, optimize = FALSE, 
                     scale = TRUE, kernel_type = "gaussian", period = NULL,
                     time_col = NULL,
                     mixed_data = FALSE, cat_columns = NULL, Xtest = NULL) {
  
  # Validate kernel_type matches available options
  kernel_type <- match.arg(kernel_type, 
                           c("gaussian", 
                             "gaussian_linear", 
                             "gaussian_quadratic",
                             "gaussian_periodic_linear", 
                             "gaussian_periodic_quadratic"))

  # Check if period is required
  needs_period <- grepl("periodic", kernel_type)
  if (needs_period && is.null(period)) {
    stop(sprintf("Error: Period parameter required for kernel type '%s'.", kernel_type))
  }
  
  # Pre-process covariates X
  X <- as.matrix(X)
  
  # Reorder columns if time_col is specified (move time column to first position)
  if (!is.null(time_col)) {
    # Check for conflict with numeric cat_columns
    if (!is.null(cat_columns) && is.numeric(cat_columns)) {
      stop("When using 'time_col', please specify 'cat_columns' by column name (character) rather than index (numeric).")
    }
    
    # Find the time column index
    if (is.character(time_col)) {
      time_idx <- which(colnames(X) == time_col)
      if (length(time_idx) == 0) {
        stop(sprintf("Specified time_col '%s' not found in X.", time_col))
      }
    } else if (is.numeric(time_col)) {
      if (time_col < 1 || time_col > ncol(X)) {
        stop(sprintf("Specified time_col index %d is out of bounds (X has %d columns).", time_col, ncol(X)))
      }
      time_idx <- time_col
    } else {
      stop("time_col must be a column name (character) or column index (numeric).")
    }
    
    # Reorder X so time column is first
    other_cols <- setdiff(1:ncol(X), time_idx)
    X <- X[, c(time_idx, other_cols), drop = FALSE]
    
    # Also reorder Xtest if provided
    if (!is.null(Xtest)) {
      Xtest <- as.matrix(Xtest)
      Xtest <- Xtest[, c(time_idx, other_cols), drop = FALSE]
    }
  }
  
  X.orig <- X
  
  # Pre-process outcome Y
  Y.init <- as.numeric(Y)
  Y.init.mean <- mean(Y.init)
  Y.init.sd <- sd(Y.init)
  Y <- scale(Y, center = Y.init.mean, scale = Y.init.sd)
  
  if (!is.null(Xtest)) {
    Xtest.orig <- as.matrix(Xtest) # Save original *unprocessed* Xtest
  } else {
    Xtest.orig <- NULL
  }
  cat_num_processed <- NULL
  cont_num_processed <- 1:ncol(X) # Default if mixed_data = FALSE
  
  ## mixed_data processing
  if(mixed_data == TRUE){
    # Xtest is passed here to get all factor levels
    X_mix <- mixed_data_processing(X.orig, cat_columns = cat_columns, Xtest = Xtest.orig) 
    X <- X_mix$X_train # use the processed X_train
    cat_num_processed <- X_mix$cat_num_processed
    cont_num_processed <- X_mix$cont_num_processed
    
    if (is.null(colnames(X))) {
      d <- ncol(X)
      colnames(X) <- paste("x", 1:d, sep = "")
    }
  } else {
    # if no mixed data, ensure X is numeric
    X <- apply(X, 2, as.numeric)
    if (is.null(colnames(X))) {
      d <- ncol(X)
      colnames(X) <- paste("x", 1:d, sep = "")
    }
  }
  X.init <- X #processed, unscaled X
  
  period_original <- period # Store original period before scaling
  
  ## Scaling & pre-kernel adjustments
  X.init.mean <- rep(0, ncol(X)) 
  names(X.init.mean) <- colnames(X)
  X.init.sd <- rep(1, ncol(X))   
  names(X.init.sd) <- colnames(X)
  
  if (scale == TRUE) { 
    # Case 1: Mixed Data (scale only continuous)
    if (!is.null(cont_num_processed) && length(cont_num_processed) > 0) {
      cont_data <- X[, cont_num_processed, drop = FALSE]
      # Calculate stats only for continuous columns
      cont_mean <- apply(cont_data, 2, mean) 
      cont_sd <- apply(cont_data, 2, sd)     
      if (sum(cont_sd == 0, na.rm=TRUE) > 0 | sum(is.na(cont_sd)) > 0) {
        stop("at least one *continuous* column in X is a constant, please remove it")
      }
      X.init.mean[cont_num_processed] <- cont_mean
      X.init.sd[cont_num_processed] <- cont_sd
      X[, cont_num_processed] <- scale(cont_data, center = cont_mean, scale = cont_sd)
      
    } else if (mixed_data == FALSE) {
      # Case 2: All Data is Continuous
      X.init.mean <- apply(X, 2, mean)
      X.init.sd <- apply(X, 2, sd)
      if (sum(X.init.sd == 0, na.rm=TRUE) > 0 | sum(is.na(X.init.sd)) > 0) {
        stop("at least one column in X is a constant, please remove the constant(s)")
      }
      X <- scale(X, center = X.init.mean, scale = X.init.sd)
    }
    
    # Scale period if provided and needed
    # Uses first continuous column (guaranteed to be time if time_col was specified)
    if (!is.null(period) && needs_period) {
      first_cont_col_idx <- NULL
      if (mixed_data == TRUE) {
        first_cont_col_idx <- cont_num_processed[1] # First from cont list
      } else {
        first_cont_col_idx <- 1 # If not mixed, it's just col 1
      }
      
      if (!is.null(first_cont_col_idx) && !is.na(first_cont_col_idx)) {
        time_col_sd <- X.init.sd[first_cont_col_idx] 
        period_scaled <- period / time_col_sd
      } else {
        warning("No continuous columns found for period scaling. Using unscaled period.")
        period_scaled <- period
      }
    } else {
      period_scaled <- NULL 
    }
    
  } else { # if scale == FALSE
    # We keep the defaults: X.init.mean = 0s and X.init.sd = 1s
    if (!is.null(period) && needs_period) {
      period_scaled <- period 
    } else {
      period_scaled <- NULL
    }
  }
  
  # apply sqrt(0.5) to categorical columns
  if(mixed_data == TRUE && !is.null(cat_num_processed)) {
    X[, cat_num_processed] <- sqrt(0.5) * X[, cat_num_processed, drop = FALSE]
  }
  
  # Optimize b by maximizing variance of K when `b=null`
  if (is.null(b)) {
    b <- getb_maxvar(X, kernel_type, period_scaled)
  }
  
  ## Calculate GP
  K <- kernel_symmetric(X, b = b, kernel_type = kernel_type, period = period_scaled)
  
  # Optimize s2, given K (with optimized b)
  if (isTRUE(optimize)) {
    opt <- gp_optimize(K = K, Y = Y)
    s2 <- opt$s2opt
  } #otherwise, user-specified s2 is given (or default s2)
  
  L <- chol(K + diag(s2, nrow(K)))
  K_inv <- Matrix::chol2inv(L)
  m <- rep(0, nrow(X)) # zero mean prior
  a <- K_inv %*% (Y - m)
  post_mean_scaled <- K %*% a
  post_cov_scaled <- K - K %*% K_inv %*% K
  prior_mean_scaled <- m
  
  # Rescale to original
  post_mean_orig <- post_mean_scaled * Y.init.sd + Y.init.mean
  post_cov_orig <- Y.init.sd^2 * post_cov_scaled
  
  results <- list(
    # Data
    X.orig = X.orig, # unprocessed original X (reordered if time_col specified)
    X.init = X.init, # processed, unscaled X
    X.init.mean = X.init.mean, # mean vector (cont only)
    X.init.sd = X.init.sd, # SD vectors (cont only)
    Y.init.mean = Y.init.mean,
    Y.init.sd = Y.init.sd,
    Y = Y,
    X = X, # final matrix X used for K (scaled for cont vars, sqrt(0.5) applied for cat vars)
    # Parameters
    b = b,
    s2 = s2,
    kernel_type = kernel_type,
    period_original = period_original,
    period_scaled = period_scaled,
    # Kernel and posterior
    K = K,
    L = L,
    alpha = as.vector(a),
    post_mean_scaled = post_mean_scaled,
    post_mean_orig = post_mean_orig,
    post_cov_scaled = post_cov_scaled,
    post_cov_orig = post_cov_orig,
    prior_mean_scaled = prior_mean_scaled,
    # Mixed data
    mixed_data = mixed_data,
    cat_columns = cat_columns, # original spec
    cat_num_processed = cat_num_processed, # processed cat indices
    cont_num_processed = cont_num_processed, # processed cont indices
    scale = scale,
    time_col = time_col, # store time_col for gp_predict
    Xcolnames = colnames(X)
  )
  return(results)
}

#' gp_predict
#'
#' to predict outcome values at testing points by feeding the results obtained from gp_train()
#'
#' @param gp a list-form object obtained from gp_train()
#' @param Xtest a data frame or a matrix of testing data set
#' @importFrom Rcpp sourceCpp
#' @return \item{Xtest_scaled}{testing data in a scaled form}
#' \item{Xtest}{the original testing data set}
#' \item{Ys_mean_scaled}{the predicted values of Y in a scaled form}
#' \item{Ys_mean_orig}{the predicted values of Y in the original scale}
#' \item{Ys_cov_scaled}{covariance of predicted Y in a scaled form}
#' \item{Ys_cov_orig}{covariance of predicted Y in the original scale}
#' \item{f_cov_orig}{covariance of target function in the original scale}
#' \item{b}{the bandwidth value obtained from gp_train()}
#' \item{s2}{the s2 value obtained from gp_train()}
#' @export

gp_predict <- function(gp, Xtest) {
  mixed_data <- gp$mixed_data
  Xtest_init <- as.matrix(Xtest) #unprocessed test data
  
  # Reorder columns if time_col was specified during training
  if (!is.null(gp$time_col)) {
    if (is.character(gp$time_col)) {
      time_idx <- which(colnames(Xtest_init) == gp$time_col)
      if (length(time_idx) == 0) {
        stop(sprintf("time_col '%s' not found in Xtest.", gp$time_col))
      }
    } else {
      time_idx <- gp$time_col
    }
    other_cols <- setdiff(1:ncol(Xtest_init), time_idx)
    Xtest_init <- Xtest_init[, c(time_idx, other_cols), drop = FALSE]
  }
  
  if(!isTRUE(mixed_data)){
    Xtest <- sweep(Xtest_init, MARGIN=2, gp$X.init.mean, FUN = "-")
    Xtest <- sweep(Xtest, MARGIN=2, gp$X.init.sd, FUN = "/")
  } else {
    Xtest_mix <- mixed_data_processing(gp$X.orig,
                                       cat_columns = gp$cat_columns,
                                       Xtest=Xtest_init)
    Xtest_processed <- Xtest_mix$X_test # one-hot encoded, unscaled test data
    
    # Check column consistency
    if (is.null(Xtest_processed)) {
      stop("Xtest processing failed, resulting matrix is NULL.")
    }
    if (ncol(Xtest_processed) != length(gp$Xcolnames)) {
      stop(sprintf("Processed Xtest has %d cols, but model was trained on %d cols. Check categorical levels.", 
                   ncol(Xtest_processed), length(gp$Xcolnames)))
    }
    
    # Reorder columns to match training order
    if (sum(colnames(Xtest_processed) != gp$Xcolnames) > 0) {
      if (!all(gp$Xcolnames %in% colnames(Xtest_processed))) {
        stop("Processed Xtest is missing columns that were present in training.")
      }
      Xtest_processed <- Xtest_processed[, gp$Xcolnames, drop = FALSE]
    }
    
    Xtest_cat <- Xtest_processed[, gp$cat_num_processed, drop = FALSE]
    Xtest_cont <- Xtest_processed[, gp$cont_num_processed, drop = FALSE]
    Xtest_cat <- sqrt(0.5) * Xtest_cat
    
    if (isTRUE(gp$scale)) {
      if (ncol(Xtest_cont) > 0) { 
        # get the correct mean/sd values from named vectors
        cont_means <- gp$X.init.mean[gp$cont_num_processed]
        cont_sds <- gp$X.init.sd[gp$cont_num_processed]
        
        Xtest_cont <- sweep(Xtest_cont, 2, cont_means, FUN = "-")
        Xtest_cont <- sweep(Xtest_cont, 2, cont_sds, FUN = "/")
      }
      # if ncol(Xtest_cont) is 0, do nothing
    }
    Xtest <- as.matrix(cbind(Xtest_cat, Xtest_cont))
    # ensure final column order matches gp$X (which is cat_cols, then cont_cols)
    colnames(Xtest) <- gp$Xcolnames
  }
  
  ## Prediction
  if (ncol(gp$X) != ncol(Xtest)) {
    stop("Internal Error: Final processed Xtest and trained X dimensions do not match.")
  }
  
  ## Calculate GP (Algorithm 2.1. in Rasmussen & Williams)
  Kss <- kernel(Xtest, Xtest, b = gp$b, kernel_type = gp$kernel_type, 
                period = gp$period_scaled) #K_{**}
  Ks <- kernel(Xtest, gp$X, b = gp$b, kernel_type = gp$kernel_type, 
               period = gp$period_scaled) #K_{*}
  
  # Compute predictive mean
  Ys_mean_scaled <- c(Ks %*% gp$alpha)
  Ys_mean_orig <- Ys_mean_scaled * gp$Y.init.sd + gp$Y.init.mean
  
  # Compute predictive covariance
  v <- solve(t(gp$L), t(Ks))
  f_cov <- Kss - crossprod(v) #cov for target function
  Ys_cov_scaled <- f_cov + diag(gp$s2, nrow(f_cov)) #cov for new observation
  # Transform back to original scale
  Ys_cov_orig <- gp$Y.init.sd^2 * Ys_cov_scaled #can be used for prediction interval
  f_cov_orig <- gp$Y.init.sd^2 * f_cov #can be used for confidence interval
  
  results <- list(
    Xtest_scaled = Xtest,
    Xtest = Xtest_init,
    #Ks = Ks, Kss = Kss,
    Ys_mean_scaled = Ys_mean_scaled,
    Ys_mean_orig = Ys_mean_orig,
    Ys_cov_scaled = Ys_cov_scaled,
    Ys_cov_orig = Ys_cov_orig,
    f_cov_orig = f_cov_orig,
    f_cov = f_cov,
    b = gp$b,
    s2 = gp$s2
  )
  
  return(results)
}


#' gp_rdd
#'
#' to perform RD using GP functions
#'
#' @param X forcing variable
#' @param Y Y vector (outcome variable)
#' @param cut cut point
#' @param alpha confidence level (default = 0.05)
#' @param b bandwidth (default = NULL)
#' @param trim a logical value indicating whether you want to do an automatic trim at a specific value of trim_k_value (default=FALSE)
#' @param trim_k_value a numerical value indicating the kernel value that you want to trim above (default = 0.1)
#' @param scale a logical value indicating whether you want to scale the covariates (default = TRUE)
#' @examples
#' n <- 100
#' tau <- 3
#' cut <- 0
#' x <- rnorm(n, 0, 1)
#' y <- rnorm(n, 0, 1) + tau*ifelse(x>cut, 1, 0)
#' gp_rdd.out <- gp_rdd(x, y, cut)
#' gp_rdd_plot(gp_rdd.out)
#' @importFrom stats sd optimize complete.cases qnorm
#' @importFrom Rcpp sourceCpp
#' @return \item{tau}{an estimated treatment effect}
#' \item{se}{the standard error of tau}
#' @export
gp_rdd <- function(X, Y, cut, alpha=0.05, b=NULL,
                   trim=FALSE, trim_k_value=0.1, scale=TRUE){
  cutpoint <- c(cut, cut)

  b_left <- b
  b_right <- b

  na.ok <- complete.cases(X) & complete.cases(Y)
  X <- as.matrix(X[na.ok])
  Y <- as.numeric(Y[na.ok])

  X_left <- X[X<cut]
  X_right <- X[X>cut]
  Y_left <- Y[X<cut]
  Y_right <- Y[X>cut]

  #trimming: use only X with kernel value <0.1
  if(isTRUE(trim) & isTRUE(scale))
    stop("If `trim==TRUE`, scale must be FALSE")
  if(isTRUE(trim)){
    #1 set b automatically using maxvarK
    b_left <- getb_maxvar(X_left)
    b_right <- getb_maxvar(X_right)

    #2 trim the sample as to remove points farther than some X value
    #whose kernel value is `trim_k_value` (default: 0.1)
    trim_at_left <- cut - sqrt(-1*b_left*log(trim_k_value))
    trim_at_right <- cut + sqrt(-1*b_right*log(trim_k_value))

    Y_left <- Y_left[X_left>trim_at_left]
    Y_right <- Y_right[X_right<trim_at_right]
    X_left <- X_left[X_left>trim_at_left]
    X_right <- X_right[X_right<trim_at_right]
  }

  ## fit GP on left side of cutoff
  gp_train_l <- gp_train(X = X_left, Y = Y_left, b=b_left, scale = scale, optimize = TRUE)
  gp_pred_l <- gp_predict(gp_train_l, cutpoint)
  ## fit GP on right side of cutoff
  gp_train_r <- gp_train(X = X_right, Y = Y_right, b=b_right, scale = scale, optimize = TRUE)
  gp_pred_r <- gp_predict(gp_train_r, cutpoint)

  ## obtain estimate, se, and CI
  pred_l <- gp_pred_l$Ys_mean_orig[1]
  pred_r <- gp_pred_r$Ys_mean_orig[1]
  tau <- pred_r - pred_l
  se <- sqrt(diag(gp_pred_l$f_cov_orig)[1] + diag(gp_pred_r$f_cov_orig)[1])
  ci <- c(lower=tau - qnorm(1-alpha/2)*se, upper=tau + qnorm(1-alpha/2)*se)

  results <- list(tau=tau, se=se, ci=ci,
                  pred_l=pred_l, pred_r=pred_r,
                  b_left=gp_train_l$b, b_right=gp_train_r$b,
                  s2_left=gp_train_l$s2, s2_right=gp_train_r$s2,
                  n_left=length(X_left), n_right=length(X_right), trim=trim,
                  X = X,
                  Y = Y,
                  gp_train_l = gp_train_l,
                  gp_train_r = gp_train_r,
                  cut = cut)

  if(trim==TRUE){
    results <- append(results,
                      c(trim_at_left=trim_at_left,
                        trim_at_right=trim_at_right))
  }
  return(results)
}

#' gp_rdd_plot
#'
#' to draw an RD plot using the results obtained from gp_rdd()
#'
#' @param gp_rdd_res a list-form results obtained from gp_rdd()
#' @param l_col a character value indicating the color of the left side of the cutoff point (default = "blue")
#' @param r_col a character value indicating the color of the right side of the cutoff point (default = "red")
#' @examples
#' library(ggplot2)
#' n <- 100
#' tau <- 3
#' cut <- 0
#' x <- rnorm(n, 0, 1)
#' y <- rnorm(n, 0, 1) + tau*ifelse(x>cut, 1, 0)
#' gp_rdd.out <- gp_rdd(x, y, cut)
#' gp_rdd_plot(gp_rdd.out) +
#'  geom_vline(xintercept = cut, lty="dashed")
#' @importFrom ggplot2 ggplot geom_point geom_line geom_ribbon theme_minimal aes
#' @return \item{gg}{an RD ggplot}
#' @export
gp_rdd_plot <- function(gp_rdd_res,
                        l_col = "blue",
                        r_col = "red"){

  ggdt <- data.frame(X = gp_rdd_res$X,
                     Y = gp_rdd_res$Y)

  Xtest_left <- seq(min(ggdt$X[ggdt$X<gp_rdd_res$cut]), gp_rdd_res$cut, length.out=100)
  gp_pred_l <- gp_predict(gp_rdd_res$gp_train_l, Xtest_left)
  l_se <- sqrt(diag(gp_pred_l$f_cov_orig))

  Xtest_right <- seq(gp_rdd_res$cut, max(ggdt$X[ggdt$X>gp_rdd_res$cut]), length.out=100)
  gp_pred_r <- gp_predict(gp_rdd_res$gp_train_r,
                          Xtest_right)
  r_se <- sqrt(diag(gp_pred_r$f_cov_orig))

  left_side <- data.frame(x = Xtest_left,
                          y = gp_pred_l$Ys_mean_orig,
                          l_lwr = gp_pred_l$Ys_mean_orig - 1.96*l_se,
                          l_upr = gp_pred_l$Ys_mean_orig + 1.96*l_se
  )

  right_side <- data.frame(x = Xtest_right,
                           y = gp_pred_r$Ys_mean_orig,
                           r_lwr = gp_pred_r$Ys_mean_orig - 1.96*r_se,
                           r_upr = gp_pred_r$Ys_mean_orig + 1.96*r_se
  )

  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(data = ggdt,
                        ggplot2::aes(X, Y), alpha = 0.5) +
    ggplot2::geom_line(data = left_side,
                       ggplot2::aes(x, y), col = l_col) +
    ggplot2::geom_ribbon(data = left_side,
                         ggplot2::aes(x, ymin = l_lwr, ymax = l_upr),
                         col = NA, alpha = 0.3, fill = l_col) +
    ggplot2::geom_line(data = right_side,
                       ggplot2::aes(x, y),
                       col = r_col) +
    ggplot2::geom_ribbon(data = right_side,
                         ggplot2::aes(x, ymin = r_lwr, ymax = r_upr),
                         col = NA, alpha = 0.3, fill = r_col) +
    ggplot2::theme_minimal()

  return(gg)
}
