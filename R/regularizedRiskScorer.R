
#' @title Fitting regularized Risk Score Models based on importance of SNPs
#' @author Mara Tews
#' @description \code{regularizedRiskScorer} is used to fit a regularized risk score models based on 
#' importance of SNPs, the number of SNPs used in the final model is calculated by an optimization
#' problem based on the models deviance less the current penalty for the number of used SNPs in the model.
#'
#' @param a [\code{integer}]\cr
#'          Number of SNPs used for inital penalty step, e.g. size of genotyping chip.
#' @param balance [\code{double}]\cr
#'          Weight of penalty.
#' @param logarithmical [\code{logical}]\cr 
#'          A logical indicating weather penalty step size should be equal (\code{FALSE}) or 
#'          logarithmic (\code{TRUE}). Default is \code{TRUE}.
#' @param data [\code{\link{data.frame}}]\cr
#'          A data frame (or object coercible by \code{\link{as.data.frame}})
#'          to a data frame) containing the variables in the model.
#' @param y.name [\code{string}]\cr
#'          A string giving the response variable name.
#' @param feature.names [\code{\link{character}}]\cr
#'          A character vector giving the feature names.
#' @param formula [\code{\link{formula}}]\cr
#'          An object of class \code{\link{formula}} (or one that can be
#'          coerced to that class): a symbolic description of the model to
#'          be fitted.
#' @param importance [\code{numeric}]\cr
#'          A vector of feature weights to be used in the fitting process.
#'          Should be a named numeric vector.
#' @param weight [\code{logical}]\cr
#'          A logical indicating weather weights should be applied
#'          (\code{TRUE}) or used to identify alleles with flipped risk
#'          allele (\code{FALSE}).
#' @param beta [\code{logical} or \code{function} or \code{number}]\cr
#'          Depends on \code{weight}: \code{TRUE}, if weights are beta
#'          coefficients; a \code{function} for transformation to beta
#'          coefficients, if \code{weight} is \code{TRUE}; a \code{number}
#'          giving a threshold as flip criteria, if \code{weight} is
#'          \code{FALSE}. See 'Details'.
#' @param ...
#'          Other parameter passed from and to other methods.
#'
#' @details Pentalty function can have equal or logarithmic step sizes. Set \code{a} as size of 
#' genotyping chip and use equal step size to determine the best amount of identical genotyping chips. 
#' Set \code{a} as a an inital step size an use logarithmic step size to regularize amount of SNPs 
#' (preferred all on one chip). 
#' 
#' @return 
#' \code{regularizedRiskScorer} returns an object of class "\code{riskScorer}". For more information see 
#' \code{\link{riskScorer}}.
#' 
#' @export
regularizedRiskScorer <- function(a, balance, logarithmical = FALSE, data, y.name, feature.names, formula, 
                                  importance, weight, beta, ...) {
  # check input and set default values if necessary
  checkmate::assertNumber(a, lower = 1, upper = length(importance))
  checkmate::assertNumber(balance, lower = 1)
  checkmate::assertLogical(logarithmical)
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
  # initialize penalty function (logarithmical or equal step size)
  penalty <- function(a, w, balance, logarithmical) {
    if(logarithmical) {
      step <- ceiling(log(((w/a)+2)/2, 1.5))
      current_penalty <- step*balance
      return(current_penalty)
    } else {
      step <- ceiling(w/a)
      current_penalty <- step*balance
      return(current_penalty)
    }
  }
  # sort importance
  imp_sorted <- importance[order(abs(importance), decreasing = TRUE)]
  # initialize optimization problem
  dev_vector <- matrix(0, nrow = length(importance), ncol = 1)
  dev_imp <- rep(0,length(importance))
  names(dev_imp) <- names(importance)
  # get intervall of possible scores
  min_resp <- sum(importance[importance<0])*2
  max_resp <- sum(importance[importance>0])*2
  interval_resp <- max_resp - min_resp
  # if necessary transform response to -1/1
  if(1 %in% data[, y.name] && 2 %in% data[, y.name]) {
    dev_truth <- 2*(as.numeric(data[, y.name])-1)-1
  } else if(0 %in% data[, y.name] && 1 %in% data[, y.name]) {
    dev_truth <- 2*(as.numeric(data[, y.name]))-1
  } else if(-1 %in% data[, y.name] && 1 %in% data[, y.name]) {
    dev_truth <- as.numeric(data[, y.name])
  } else {
    stop("Status of disease should be given as -1/1, 0/1 or 1/2.")
  }
  dev_data <- as.matrix(data[, names(importance)])
  dev_modelling <- function(w) {
    # only use importance 1:w, other values = 0
    dev_imp[names(imp_sorted[1:w])] <- imp_sorted[1:w]
    # return probability for class 1 
    resp_temp <- dev_data%*%dev_imp
    resp_temp <- (resp_temp-min_resp)/interval_resp 
    # get devianz for each model
    dev_vector[w] <- sum(-log(1+exp(-2*dev_truth*(2*resp_temp-1))))
  }
  mod_deviance <- sapply(
    X = 1:length(importance), 
    FUN = dev_modelling
    )
  # get penalty
  mod_penalty <- sapply(
    X = 1: length(mod_deviance),
    FUN = function (s) penalty(
      w = s, 
      a = a, 
      balance = balance, 
      logarithmical = logarithmical
      )
    )
  # solve optimization problem
  best_model <- which.max(mod_deviance - mod_penalty)
  # get final model for prediction
  fin_imp <- imp_sorted[1:best_model]
  fin_feature.names <- names(fin_imp) 
  fin_y.name <- y.name
  fin_data <- dplyr::select(data, fin_feature.names, y.name)
  fin_model <- riskScorer(
    data = fin_data, 
    y.name = fin_y.name, 
    feature.names = names(fin_imp),
    importance = fin_imp,
    weight = weight, 
    beta = TRUE
    )
  fin_model$balance <- balance
  class(fin_model) <- c(class(fin_model), "regularizedRiskScorer")
  return(fin_model)
}
