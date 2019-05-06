library(mlr)

# definition of the learner

makeRLearner.classif.regularizedRiskScorer = function() {
  makeRLearnerClassif(
    cl = "classif.regularizedRiskScorer",
    package = "riskScoreR",
    par.set = makeParamSet(
      makeNumericLearnerParam(id = "a", lower = 1, upper = Inf, when = "train", tunable = TRUE),
      makeNumericLearnerParam(id = "balance", lower = 1, upper = Inf, when = "train", tunable = TRUE),
      makeLogicalLearnerParam(id = "weight", when = "train", tunable = FALSE),
      makeNumericLearnerParam(id = "beta", lower = 0, upper = Inf, when = "train", requires = (weight == TRUE), tunable = FALSE),
      makeLogicalLearnerParam(id = "beta", when = "train", requires = (weight == TRUE), tunable = FALSE),
      makeFunctionLearnerParam(id = "beta", when = "train", requires = (weight == FALSE), tunable = FALSE),
      makeDiscreteLearnerParam(id = "type", default = "score", values = c("score", "prob", "response"), when = "both", tunable = FALSE)),
    properties = c("twoclass", "numerics", "prob"), 
    name = "Regularized Risk Scorer",
    short.name = "regRS")
  }

# creating the training function of the learner

trainLearner.classif.regularizedRiskScorer = function(.learner, .task, .subset, .weights = NULL, ...) {
  riskScoreR::regularizedRiskScorer(data = getTaskData(.task, .subset),
                                    importance = .weights,
                                    y.name = getTaskTargetNames(.task),
                                    feature.names = getTaskFeatureNames(.task),
                                    formula = getTaskFormula(.task))
  }

# creating the prediction method

predictLearner.classif.regularizedRiskScorer = function(.learner, .model, .newdata, type, ...) {
  p = predict(risk.scorer = .model$learner.model, 
              newdata = .newdata, 
              type = type)
  }




