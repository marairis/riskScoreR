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
#' @param weights [\code{numeric}]\cr
#'                A vector of weights to be used in the fitting process. Should
#'                be a named numeric vector.
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
#'
#' @details Generally, weights are assumed to be beta coefficients from a GWAS
#' analysis. In this case, if \code{weight} is \code{TRUE}, the number of risk
#' alleles are multiplied with it's beta coefficients to get a risk score. If no
#' beta coefficients are available, a function to transform the given weights to
#' beta coefficient equivalents is needed to be specified by \code{beta} (might
#' be \code{\link{identity}}). If \code{weight} is \code{FALSE} no weighting of
#' risk alleles is performed. In this case a threshold is needed to find flipped
#' risk alleles by the given weights. If beta coefficients are given (default)
#' this threshold is set to 0. It would be 1 if the weights are odds ratios.
#' The parameter \code{beta} can be set to any finite number.
#'
#' @return \code{riskScorer} returns an object of class "\code{riskScorer}". An
#' object of class "\code{riskScorer}" is a list containing the following
#' components:
#' \describe{
#'  \item{\code{riskModel}}{An object of class \code{\link{data.table}}
#'                          specifying the risk model. Contains a column
#'                          \code{variable} with the name of the variant, a
#'                          column \code{weight} specifying its weight and an
#'                          optional column \code{flip} indicating if the risk
#'                          alleles of this variant have to be flipped.}
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
riskScorer <- function(formula, data, weights, weight, beta = TRUE) {

  # Check input
  checkmate::assertClass(f <- as.formula(formula), "formula")
  checkmate::assertDataFrame(data)
  checkmate::assertNumeric(weights, any.missing = FALSE)
  checkmate::assertNamed(weights)
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

  if(any(as.character(f[[3]]) != ".")) {
    data <- data.table::as.data.table(model.frame(f, data))
  }
  data.table::setnames(data, colnames(data), make.names(colnames(data)))

  riskModel <- weights[intersect(names(weights), colnames(data))]

  riskModel <- data.table::data.table(variable = names(riskModel),
                                      weight = riskModel,
                                      key = "variable")

  if(weight) {
    riskModel[, weight := beta(weight)]
  } else {
    # find variants to flip risk allele
    riskModel[, flip := weight < beta]
    riskModel[, weight := 1]
  }

  out <- list()
  class(out) <- "riskScorer"
  out$riskModel <- riskModel
  out$beta <- beta
  out$weight <- weight
  out$formula <- f

  return(out)
}

#' @title Predict Method for riskScoreR Fits
#'
#' @description Obtains predictions from a fitted risk score model object.
#'
#' @param riskScorer  [\code{\link{riskScorer}}]\cr
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
#'
#' @details Predicted probabilities are obtained by rescaling the distribution
#' of scores to \eqn{[0,1]}. Predicted class labels are obtained by rounding the
#' predicted probabilites.
#'
#' @return An object of class "\code{data.frame}" is returned with two columns.
#' First column is the id. If \code{type = "score"} (un)weighted risk scores; if
#' \code{type = "respons"} predicted class labels; if \code{type = "prob"}
#' predicted class probabilities are returend in the second column.
#'
#' @import data.table
#'
#' @export
predict.riskScorer <- function(riskScorer, newdata, type = "score") {

  # Check input
  checkmate::assertClass(riskScorer, "riskScorer")
  if(!checkmate::testClass(newdata, "data.table")) {
    checkmate::assertDataFrame(newdata)
    rn <- rownames(newdata)
    newdata <- data.table::as.data.table(newdata)
    rownames(newdata) <- rn
  }
  if(any(as.character(riskScorer$formula[[3]]) != ".")) {
    newdata <- data.table::as.data.table(model.frame(riskScorer, newdata))
  }
  if(!"id" %in% colnames(newdata)) {
    data.table::set(newdata, j="id", value=rownames(newdata))
  }
  checkmate::assertChoice(type, c("score", "response", "prob"))

  target <- as.character(riskScorer$formula[[2]])
  target.levels <- levels(factor(newdata[[target]]))
  if(nlevels(newdata[[target]]) != 2) {
    stop("Only two class problems are allowed!")
  }

  newdata <- melt(newdata,
                  id.vars = c("id", target))
  data.table::setkey(newdata, variable)

  newdata <- riskScorer$riskModel[newdata]

  if(!riskScorer$weight){
    newdata[flip == TRUE, value := 2-value]
  }

  scores <- newdata[, .(score=sum(value*weight)), by = id]

  if(type == "score") {
    return(as.data.frame(scores))
  } else {
    prob <- scores[, .(id = id,
                       good = 1-scales::rescale(score),
                       bad = scales::rescale(score))]
    if(type == "response") {
      return(as.data.frame(prob[, .(id = id, response=factor(target.levels[1+round(prob[["bad"]])]))]))
    } else {
      data.table::setnames(prob, c("bad", "good"), c(target.levels[[2]], target.levels[[1]]))
      return(as.data.frame(prob))
    }
  }

}