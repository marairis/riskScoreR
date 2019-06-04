#' @title Feature selection filter
#' 
#' @description Selects optimal number of features to fit a riskScorer model.
#' 
#' @import mlr
#' @export
makeFilter(
  name = "feature_selection.filter",
  desc = "Selects a subset of relevant features for use in model construction. ",
  pkg = "riskScoreR",
  supported.tasks = "classif",
  supported.features = "numerics",
  fun = function(task, nselect, learner, decreasing = TRUE, importance, ...){
    return(sort(importance, decreasing = decreasing))
    }
  )