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

#' print method for objects of class gpss
#'
#' @param x gpss object
#' @examples
#' library(gpssinternal)
#' data(lalonde)
#' # categorical variables must be encoded as factors
#' dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))
#' # train and test sets
#' idx <- sample(seq_len(nrow(dat)), 500)
#' dat_train <- dat[idx, ]
#' dat_test <- dat[-idx, ]
#' # sample of data for speed
#' mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat_train)
#' print(mod)
#' @export
print.gpss <- function(x, ...) {
  cat("gpss object\n")
  cat("===========\n")
  cat("Formula: "); print(attr(x, "formula"))
  cat("Kernel type:", x$kernel_type, "\n")
  cat("Bandwidth (b):", x$b, "\n")
  cat("Noise variance (s2):", x$s2, "\n")
  
  # Periodic kernel info
  if (!is.null(x$period_original)) {
    cat("Period (original):", x$period_original, "\n")
    cat("Period (scaled):", x$period_scaled, "\n")
    if (!is.null(x$time_col)) {
      cat("Time column:", x$time_col, "\n")
    }
  }
  
  # Mixed data info
  cat("Mixed data:", x$mixed_data, "\n")
  if (isTRUE(x$mixed_data)) {
    cat("Categorical columns:", x$cat_columns, "\n")
  }
  
  cat("\nUse `summary()` for more details.\n")
  invisible(x)
}

#' summary method for objects of class gpss
#'
#' @param object gpss object
#' @examples
#' library(gpssinternal)
#' data(lalonde)
#' # categorical variables must be encoded as factors
#' dat <- transform(lalonde, race_ethnicity = factor(race_ethnicity))
#' # train and test sets
#' idx <- sample(seq_len(nrow(dat)), 500)
#' dat_train <- dat[idx, ]
#' dat_test <- dat[-idx, ]
#' # sample of data for speed
#' mod <- gpss(re78 ~ nsw + age + educ + race_ethnicity, data = dat_train)
#' summary(mod)
#' @export
summary.gpss <- function(object, ...) {
  cat("GPSS Model Summary\n")
  cat("==================\n\n")
  
  # Basic info
  cat("Model Information\n")
  cat("-----------------\n")
  cat("Formula: "); print(attr(object, "formula"))
  cat("Observations:", nrow(object$Y), "\n")
  cat("Covariates:", ncol(object$X), "\n")
  cat("Scaled:", object$scale, "\n\n")
  
  # Kernel info
  cat("Kernel Settings\n")
  cat("---------------\n")
  cat("Type:", object$kernel_type, "\n")
  cat("Bandwidth (b):", object$b, "\n")
  cat("Noise variance (s2):", object$s2, "\n")
  
  if (!is.null(object$period_original)) {
    cat("Period (original):", object$period_original, "\n")
    cat("Period (scaled):", object$period_scaled, "\n")
    if (!is.null(object$time_col)) {
      cat("Time column:", object$time_col, "\n")
    } else {
      cat("Time column: first continuous column (default)\n")
    }
  }
  cat("\n")
  
  # Mixed data info
  cat("Data Composition\n")
  cat("----------------\n")
  cat("Mixed data:", object$mixed_data, "\n")
  if (isTRUE(object$mixed_data)) {
    cat("Categorical columns:", paste(object$cat_columns, collapse = ", "), "\n")
    cat("Processed categorical indices:", paste(object$cat_num_processed, collapse = ", "), "\n")
    cat("Processed continuous indices:", paste(object$cont_num_processed, collapse = ", "), "\n")
  }
  cat("\n")
  
  # Usage hints
  cat("Usage\n")
  cat("-----\n")
  cat("Fitted values: fit$post_mean_orig\n")
  cat("Standard errors: sqrt(diag(fit$post_cov_orig))\n")
  cat("Predictions: predict(fit, newdata)\n")
  
  invisible(object)
}
