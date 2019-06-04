#' @title Multi-criteria tuning for regularizedRiskScorer
#'
#' @description Multi-criteria tuning for balance and number of used features for regularizedRiskScorer model. 
#'
#' @param a [\code{integer}]\cr
#'          Number of SNPs used for inital penalty step, e.g. size of genotyping chip.
#' @param min_balance [\code{numeric}]\cr
#'          Lower limit for balance.
#' @param max_balance [\code{numeric}]\cr
#'          Upper limit for balance.
#' @param k [\code{numeric}]\cr
#'          Number of cross-validation steps. Default is 10.
#' @param tune_strategy [\code{string}]\cr
#'          Pass mlr \code{Tune Controll} object to spezify search strategy and resolution.
#'          Search strategy used in cross-validation to optimize \code{balance}. One of those implemented in
#'          \code{\link{mlr}} package. Default is grid search with \code{5L} resolution.#' 
#' @param logarithmical [\code{logical}]\cr 
#'          A logical indicating weather penalty step size should be equal (\code{FALSE}) or logarithmic (\code{TRUE}).
#'          Default is \code{TRUE}.
#' @param data [\code{\link{data.frame}}]\cr
#'          A data frame (or object coercible by \code{\link{as.data.frame}})
#'          to a data frame) containing the variables in the model.
#' @param feature.names [\code{\link{character}}]\cr
#'          A character vector giving the feature names.
#' @param y.name [\code{string}]\cr
#'          A string giving the response variable name.
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
#'          \code{FALSE}.
#'  
#' @param ...
#'          Other parameter passed from and to other methods.
#'
#' @export
multicrit.regularizedRiskScorer <- function(a, min_balance, max_balance, k, tune_strategy, logarithmical, 
                                      data, feature.names, y.name, formula, importance, weight, beta, ...) {
  # check input and set default values if necessary
  checkmate::assertNumber(a, lower = 1, upper = length(importance))
  checkmate::assertNumeric(min_balance, lower = 1, upper = max_balance)
  checkmate::assertNumeric(max_balance, lower = min_balance)
  checkmate::assertNumber(k, lower = 1)
  if(missing(tune_strategy)) {
    tune_strategy <- mlr::makeTuneMultiCritControlGrid(resolution = 5L)
  } else {
    checkmate::assertClass(tune_strategy, classes = "TuneMultiCritControl")
  }
  if(missing(logarithmical)) {
    checkmate::assertLogical(logarithmical <- FALSE)
  } else {
    checkmate::assertLogical(logarithmical)
  }
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
  # initialize tuning
  rrs_training_task <- mlr::makeClassifTask(
    id = "training", 
    data = data,
    target = y.name
    )
  rrs_learner <- mlr::makeLearner(
    "classif.regularizedRiskScorer",
    predict.type = "prob", 
    a = a,
    importance = importance,
    logarithmical = logarithmical,
    weight = weight,
    beta = beta
    )
  rrs_param <- ParamHelpers::makeParamSet(
    makeNumericParam(
      "balance", 
      lower = min_balance, 
      upper = max_balance
      )
    )
  rrs_desc <- mlr::makeResampleDesc(
    method = "CV",
    predict = "test",
    iters = k
    )
  rrs_tune_control <- tune_strategy
  # tune balance
  rrs_tuned <- mlr::tuneParamsMultiCrit(
    learner = rrs_learner,
    task = rrs_training_task,
    resampling = rrs_desc,
    measures = list(auc, my.bal), 
    par.set = rrs_param,
    control = rrs_tune_control
    )
  }