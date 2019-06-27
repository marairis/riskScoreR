#' @title Feature selection filter
#' @author Mara Tews
#' @description Creates a filter to select the optimal number of features to fit a riskScorer model.
#' 
#' @import mlr
#' 
#' @export
featureSelection.filter <- makeFilter(
  name = "featureSelection.filter",
  desc = "Selects a subset of relevant features for use in model construction.",
  pkg = "riskScoreR",
  supported.tasks = "classif",
  supported.features = "numerics",
  fun = function(task, nselect, learner, decreasing = TRUE, importance, ...){
    return(sort(importance, decreasing = decreasing))
    }
)
