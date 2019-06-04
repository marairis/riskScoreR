#' @import mlr
#' @import ParamHelpers
#' @import stats
#--- riskScorer ---
# definition of the learner
#' @export
makeRLearner.classif.riskScorer <- function() {
  makeRLearnerClassif(
    cl = "classif.riskScorer",
    package = "riskScoreR",
    par.set = makeParamSet(
      makeLogicalLearnerParam(id = "weight", when = "train", tunable = TRUE),
      makeUntypedLearnerParam(id = "beta", when = "both", tunable = FALSE),
      makeNumericVectorLearnerParam(id = "importance", tunable = FALSE)
      ),
    properties = c("twoclass", "numerics", "prob"),
    name = "Risk Scorer",
    short.name = "rs"
    )
  }
# creating the training function of the learner
#' @export
trainLearner.classif.riskScorer <- function(.learner, .task, .subset, .weights = NULL, ...) {
  riskScoreR::riskScorer(
    data = getTaskData(.task, .subset),
    formula = getTaskFormula(.task), 
    ...
    )
  }
# creating the prediction method
#' @export
predictLearner.classif.riskScorer <- function(.learner, .model, .newdata, ...) {
  # fallunterscheidung nach type
  p = predict(
    object = .model$learner.model, 
    newdata = .newdata, 
    type = .learner$predict.type
    )
  if (.learner$predict.type == "response") {
    p[[1]]
  } else {
    as.matrix(p)
  }
  }
# register the new S3 methodes
registerS3method("makeLearner", "riskScorer", makeRLearner.classif.riskScorer)
registerS3method("trainLearner", "riskScorer", trainLearner.classif.riskScorer)
registerS3method("predictLearner", "riskScorer", predictLearner.classif.riskScorer)

#--- regularizedRiskScorer ---
# definition of the learner
#' @export
makeRLearner.classif.regularizedRiskScorer <- function() {
  makeRLearnerClassif(
    cl = "classif.regularizedRiskScorer",
    package = "riskScoreR",
    par.set = makeParamSet(
      makeNumericLearnerParam(id = "a", lower = 1, upper = Inf, when = "train", tunable = FALSE),
      makeNumericLearnerParam(id = "balance", lower = 1, upper = Inf, when = "train", tunable = TRUE),
      makeLogicalLearnerParam(id = "weight", when = "train", tunable = TRUE),
      makeLogicalLearnerParam(id = "logarithmical", when = "train", tunable = FALSE),
      makeUntypedLearnerParam(id = "beta", when = "both", tunable = FALSE),
      makeNumericVectorLearnerParam(id = "importance", tunable = FALSE)
      ),
    properties = c("twoclass", "numerics", "prob"), 
    name = "Regularized Risk Scorer",
    short.name = "rrs"
    )
  }
# creating the training function of the learner
trainLearner.classif.regularizedRiskScorer <- function(.learner, .task, .subset, .weights = NULL, ...) {
  riskScoreR::regularizedRiskScorer(
    data = getTaskData(.task, .subset),
    formula = getTaskFormula(.task), 
    ...
    )
  }
# creating the prediction method
predictLearner.classif.regularizedRiskScorer <- function(.learner, .model, .newdata, ...) {
  # fallunterscheidung nach type
  p = predict(
    object = .model$learner.model, 
    newdata = .newdata, 
    type = .learner$predict.type
    )
  if (.learner$predict.type == "response") {
    p[[1]]
  } else {
    as.matrix(p)
  }
  }
# register the new S3 methodes
registerS3method("makeLearner", "regularizedRiskScorer", makeRLearner.classif.regularizedRiskScorer)
registerS3method("trainLearner", "regularizedRiskScorer", trainLearner.classif.regularizedRiskScorer)
registerS3method("predictLearner", "regularizedRiskScorer", predictLearner.classif.regularizedRiskScorer)

#--- balance measure ---
# calculate balance
my.bal.fun <- function(task, model, pred, feats, extra.args) {
  return(model$learner.model$balance)
}
# generate measure object
my.bal <- makeMeasure(
  id = "my.bal",
  name = "My Balance Measure",
  properties = c("classif"),
  minimize = FALSE,
  best = Inf,
  worst = 0,
  fun = my.bal.fun
  )
