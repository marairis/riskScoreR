
#' @title Optimization of number of used SNPs.
#'
#' @description Optimizes number of used SNPs w used in final model based on penaly function.
#'
#' @param a [\code{integer}]\cr
#'          Initial number of SNPs.
#' @param data [\code{\link{data.frame}}]\cr
#'          A data frame (or object coercible by \code{\link{as.data.frame}})
#'          to a data frame) containing the variables in the model.
#' @param grid.res [\code{integer}]\cr
#'          Resolution of the grid for w, e.g. number of equal-sized steps. Default is 10.
#' @param importance [\code{numeric}]\cr
#'          A vector of feature weights to be used in the fitting process.
#'          Should be a named numeric vector.
#' @param k [\code{integer}]\cr
#'          Number of folds used in cross-validation to find best w. Default is 10.
#' @param balance [\code{double}]\cr
#'          Weight of penalty.
#' @param ...
#'          Other parameter passed from and to other methods.
#'
#' @return Optimized number of SNPs w to use in final model.
#'
#' @export
#'
regularizedRiskScorer <- function(a, balance, data, importance, k, grid.res, weight, formula, beta, ...) {
  # (i) check input and set default values
  checkmate::assertNumber(a)
  checkmate::assertNumber(balance)
  checkmate::assertDataFrame(data, col.names = "named")
  checkmate::assertNumeric(importance, any.missing = FALSE)
  checkmate::assertNamed(importance)

  if(!missing(formula)) {
    checkmate::assertClass(f <- stats::as.formula(formula), "formula")
    y.name <- all.vars(f[-3])
    feature.names <- all.vars(f[-2])
    if("." %in% feature.names) {
      feature.names <- setdiff(colnames(data), y.name)
    }
  } else {
    if(missing(y.name) && missing(feature.names)) {
      stop("If not providing a formula, 'y.name' and 'feature.names' must be given!")
    }
  }

  # sort importance
  imp_sorted <- importance[order(abs(importance), decreasing = TRUE)]

  # initialize penalty function
  penalty <- function(a, w, balance) {
    step <- ceiling(log(((w/a)+2)/2, 1.5))
    current_penalty <- step*balance
    return(current_penalty)
  }

  # (ii) initilize optimization problem (min)
  dev_vector <- matrix(0, nrow = length(importance), ncol = 1)
  dev_imp <- rep(0,length(importance))
  names(dev_imp) <- names(importance)
  # get intervall of possible scores
  min_resp <- sum(importance[importance<0])*2
  max_resp <- sum(importance[importance>0])*2
  interval_resp <- max_resp - min_resp
  # transform 0/1 to -1/1
  dev_truth <- (2*data[, y.name])-1
  dev_data <- as.matrix(data[, names(importance)])

  dev_modelling <- function(w) {
    # only use importance 1:w, other values = 0
    dev_imp[names(imp_sorted[1:w])] <- imp_sorted[1:w]
    # return probability
    resp_temp <- dev_data%*%dev_imp
    resp_temp <- resp_temp/interval_resp
    # get devianz for each model
    dev_vector[w] <- sum(-log(1+exp(-2*dev_truth*(2*resp_temp-1))))
  }
  # returns vector with deviance for each model
  mod_deviance <- sapply(X = 1:length(importance),
                       FUN = dev_modelling)

  # get penalty
  mod_penalty <- sapply(X = 1: length(mod_deviance),
                        FUN = function (s) penalty(w = s, a = a, balance = balance))

  # solve optimization problem: max(devianz - penalty)
  best_model <- which.max(mod_deviance - mod_penalty)

  # get final model for prediction
  fin_imp <- imp_sorted[1:best_model]
  fin_feature.names <- names(fin_imp)
  fin_y.name <- y.name
  fin_data <- dplyr::select(data, fin_feature.names, y.name)

  fin_model <- riskScorer(data = fin_data,
                          y.name = fin_y.name,
                          feature.names = names(fin_imp),
                          importance = fin_imp,
                          weight = weight,
                          beta = TRUE)
  return(fin_model)
}