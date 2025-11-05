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
print.gpss <- function(x) {
  cat("gpss object contains information, including...\n")
  cat("formula: "); print(attr(x, "formula"))
  cat("b (bandwidth):", x$b, "\n")
  cat("s2 (noise variance):", x$s2, "\n")
  cat("mixed data (containing a categorical variable?):", x$mixed_data, "\n")
  if(x$mixed_data == TRUE){
    cat("categorical columns:", x$cat_columns,"\n")
    cat("categorical column numbers:", x$cat_num,"\n\n")
  }else{
    cat("\n")
  }
  cat("Other available information: posterior mean (scaled and original), posterior covariance (scaled and original), a kernel matrix of X, original values of Y and X, etc.", "\n\n")
  cat("For more detailed information, please use `summary(gpss object)`")
}


