
flip <- NULL
response <- NULL

#' @title Fitting Risk Score Models
#'
#' @description \code{riskScorer} is used to fit a risk score models, specified
#' by giving a symbolic description of the linear predictor and weights from a,
#' for example GWAS, analysis.
#'
#' @param formula [\code{\link{formula}}]\cr
#'                An object of class \code{\link{formula}} (or one that can be
#'                coerced to that class): a symbolic description of the model to
#'                be fitted.
#' @param data    [\code{\link{data.frame}}]\cr
#'                A data frame (or object coercible by \code{\link{as.data.frame}})
#'                to a data frame) containing the variables in the model.
#' @param y.name  [\code{string}]\cr
#'                A string giving the response variable name.
#' @param feature.names [\code{\link{character}}]\cr
#'                      A character vector giving the feature names.
#' @param importance [\code{numeric}]\cr
#'                A vector of feature weights to be used in the fitting process.
#'                Should be a named numeric vector.
#' @param weight  [\code{logical}]\cr
#'                A logical indicating weather weights should be applied
#'                (\code{TRUE}) or used to identify alleles with flipped risk
#'                allele (\code{FALSE}).
#' @param beta    [\code{logical} or \code{function} or \code{number}]\cr
#'                Depends on \code{weight}: \code{TRUE}, if weights are beta
#'                coefficients; a \code{function} for transformation to beta
#'                coefficients, if \code{weight} is \code{TRUE}; a \code{number}
#'                giving a threshold as flip criteria, if \code{weight} is
#'                \code{FALSE}. See 'Details'.
#' @param ...     Other parameter passed from and to other methods.
#'
#' @details Generally, importance values/feature weights are assumed to be beta
#' coefficients from a GWAS analysis. In this case, if \code{weight} is
#' \code{TRUE}, the number of risk alleles are multiplied with it's beta
#' coefficients to get a risk score. If no beta coefficients are available, a
#' function to transform the given weights to beta coefficient equivalents is
#' needed to be specified by \code{beta} (might be \code{\link{identity}}). If
#' \code{weight} is \code{FALSE} no weighting of risk alleles is performed. In
#' this case a threshold is needed to find flipped risk alleles by the given
#' weights. If beta coefficients are given (default) this threshold is set to 0.
#' It would be 1 if the weights are odds ratios. The parameter \code{beta} can
#' be set to any finite number.
#'
#' @return \code{riskScorer} returns an object of class "\code{riskScorer}". An
#' object of class "\code{riskScorer}" is a list containing the following
#' components:
#' \describe{
#'  \item{\code{risk_model}}{An object of class \code{\link{data.table}}
#'                          specifying the risk model. Contains a column
#'                          \code{variable} with the name of the variant, a
#'                          column \code{weight} specifying its weight and a
#'                          column \code{flip} indicating if the risk alleles of
#'                          this variant have to be flipped.}
#'  \item{\code{beta}}{Depends on \code{weight}: \code{TRUE}, if weights are beta
#'                     coefficients; a \code{function} for transformation to
#'                     beta coefficients, if \code{weight} is \code{TRUE}; a
#'                     \code{number} giving a threshold as flip criteria, if
#'                     \code{weight} is \code{FALSE}.}
#'  \item{\code{weight}}{A logical indicating if weighting of risk alleles has
#'                       to be done.}
#'  \item{\code{formula}}{An object of class \code{\link{formula}}: a symbolic
#'                        description of the model to be fitted.}
#' }
#'
#' @import data.table
#'
#' @export
riskScorer <- function(formula, data, y.name, feature.names, importance, weight, beta = TRUE, ...) {

  # Check input
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
  checkmate::assertNumeric(importance, any.missing = FALSE)
  checkmate::assertNamed(importance)
  checkmate::assertDataFrame(data, col.names = "named")
  if(!missing(y.name)) {
    checkmate::assertString(y.name, na.ok = FALSE, null.ok = FALSE)
    checkmate::assertSubset(y.name, colnames(data), empty.ok = FALSE)
    y <- factor(data[[y.name]])
    checkmate::assertFactor(y, n.levels = 2)
    target <- y.name
    target_levels <- levels(y)
  }
  if(!missing(feature.names)) {
    checkmate::assertCharacter(feature.names, all.missing = FALSE, any.missing = FALSE, null.ok = FALSE)
    checkmate::assertSubset(feature.names, names(importance), empty.ok = FALSE)
  }
  checkmate::assertLogical(weight)
  if(weight) {
    if(checkmate::testLogical(beta)) {
      if(beta) {
        beta <- base::identity
      } else {
        stop("If weights are not beta coefficients and weight is TRUE, you have to specify a function for transformation!")
      }
    } else {
      checkmate::assertFunction(beta)
    }
  } else {
    if(checkmate::testLogical(beta)) {
      if(beta) {
        beta <- 0
      } else {
        stop("If weights are not beta coefficients and weight is FALSE, you have to specify a threshold as flip criteria!")
      }
    } else {
      checkmate::assertNumber(beta, finite = TRUE)
    }
  }

  risk_model <- importance[intersect(names(importance), feature.names)]

  risk_model <- data.table::data.table(variable = names(risk_model),
                                      weight = risk_model,
                                      key = "variable")

  if(weight) {
    risk_model[, flip := FALSE]
    risk_model[, weight := beta(weight)]
  } else {
    # find variants to flip risk allele
    risk_model[, flip := weight < beta]
    risk_model[, weight := 1]
  }

  out <- list()
  class(out) <- "riskScorer"
  out$call <- match.call()
  out$risk_model <- risk_model
  out$beta <- beta
  out$weight <- weight
  out$target <- target
  out$target_levels <- target_levels

  return(out)
}

#' @title Predict Method for riskScoreR Fits
#'
#' @description Obtains predictions from a fitted risk score model object.
#'
#' @param risk.scorer [\code{\link{riskScorer}}]\cr
#'                    A fitted object of class "\code{riskScorer}".
#' @param newdata     [\code{data.frame}]\cr
#'                    A data frame in which to look for variables with which to
#'                    predict.
#' @param type        [\code{string}]\cr
#'                    The type of prediction required. The default is the
#'                    aggregated risk score, i.e. the (un)weighted sum of risk
#'                    alleles. The alternative "\code{response}" gives the
#'                    predicted binary class labels. The alternative
#'                    "\code{prob}" gives the  predicted probabilities of two
#'                    classes.
#' @param ...         Other parameter passed from and to other methods.
#'
#' @details Predicted probabilities are obtained by rescaling the distribution
#' of scores to \eqn{[0,1]}. Predicted class labels are obtained by rounding the
#' predicted probabilites.
#'
#' @return Depends on type: if \code{type = "score"} a \code{data.frame} with
#' one column holding the raw scores, i.e. (un)weighted sum of risk alleles; if
#' \code{type = "prob"} a \code{data.frame} with two columns, named as the two
#' factor levels in the response variable of the \code{riskScorer} object, giving
#' the probabilities of the respective factor level; if \code{type = "response"}
#' a \code{data.frame} holding the factor level with highest probability.
#'
#' @import data.table
#'
#' @export
predict.riskScorer <- function(risk.scorer, newdata, type = "score", ...) {

  # Check input
  checkmate::assertClass(risk.scorer, "riskScorer")
  checkmate::assertDataFrame(newdata)
  checkmate::assertChoice(type, c("score", "response", "prob"))

  # Select variables in risk model
  if(!all(risk.scorer$risk_model$variable %in% colnames(newdata))) {
    stop("Some variables listed in the given risk model are not present in newdata!")
  } else {
    newdata <- as.matrix(subset(newdata, select = risk.scorer$risk_model$variable))
  }

  if(risk.scorer$weight) {
    # don't flip genotypes
    scores <- newdata %*% risk.scorer$risk_model$weight
  } else {
    # flip genotypes
    scores <- abs(((rep_len(2, length.out = nrow(newdata))) %*% t(risk.scorer$risk_model$flip)) - newdata) %*% risk.scorer$risk_model$weight
  }

  if(type == "score") {
    return(data.frame(score = scores))
  } else {
    # High score is bad, low score is good
    prob <- data.frame(good = 1-scales::rescale(scores),
                       bad = scales::rescale(scores))
    if(type == "response") {
      prob$response <- factor(risk.scorer$target_levels[1+round(prob$bad)])
      return(subset(prob, select = response))
    } else {
      colnames(prob) <- risk.scorer$target_levels
      return(as.data.frame(prob))
    }
  }

}