
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

optimized_measure <- function(a, balance, data, importance, k, grid.res, weight, formula, beta, ...) {
  ## die anderen variablen auch noch checken?
  # (i) check input and set default values
  checkmate::assertNumber(a)
  checkmate::assertNumber(balance)
  checkmate::assertDataFrame(data, col.names = "named")
  checkmate::assertNamed(importance)
  if(missing(k)) {
    checkmate::assertNumber(k <- 10)
  }
  if(missing(grid.res)) {
    checkmate::assertNumber(grid.res <- 10)
  }

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

  # initialize penalty function
  penalty <- function(a, w, balance, data) {
    step <- ceiling(log(((w/a)+2)/2, 1.5))
    current_penalty <- step*balance
    # transform into [0,1] based on maximum number of w
    max_penalty <- (ceiling(log((((ncol(data)-1)/a)+2)/2, 1.5)))*balance
    penalty2auc <- current_penalty/max_penalty
    return(penalty2auc)
  }

  # (ii) get ready for cross-validation
  n_data <- nrow(data)
  # ranking snps based on impotance (from gwas)
  variables_ranked <- sort(importance, decreasing = TRUE)
  # shuffle data for crossvalidation
  perm_idx <- sample(nrow(data))
  # calculate possible values for w, ncol(data)-1 wegen state
  w_candidates <- sapply(list(1:grid.res), function(s) ((ncol(data)-1)/grid.res)*s)
  # describe modelling function: does a cross-validation and returns performance for current w

  modelling <- function(j) {
    cv_auc_sum <- 0
    for (i in 0:k-1) {
      # select train and test data
      cv_test <- data[perm_idx[((n_data/k)*i+1):((n_data/k)*(i+1))], ]
      cv_train <- data[-perm_idx[((n_data/k)*i+1):((n_data/k)*(i+1))], ]

      # set up model
      # create parameters adjusted to w
      ## noch allgemein formulieren
      #cv_formula <- as.formula(paste("state ~", paste(paste0("x", 1:w_candidates[j]), collapse = "+")))
      cv_y.name <- y.name
      cv_feature.names <- feature.names[1:w_candidates[j]]
      cv_importance <- unlist(importance_data)[1:w_candidates[j]]
      cv_data <- cv_train[, 1:w_candidates[j]]
      cv_data$state <- cv_train$state

      cv_model <- riskScorer(data = cv_data,
                             y.name = cv_y.name,
                             feature.names = cv_feature.names,
                             importance = cv_importance,
                             weight,
                             beta = TRUE)

      # validate on test data
      cv_newdata <- cv_test[,1:w_candidates[j]]
      colnames(cv_newdata) <- cv_feature.names
      cv_newdata$state <- cv_test$state
      cv_prediction <- predict.riskScorer(risk.scorer = cv_model,
                                          newdata = cv_newdata,
                                          type = "response")

      # get auc
      cv_auc <- mlr::measureAUC(probabilities = cv_prediction,
                                truth = cv_test$state,
                                negative = 0,
                                positive = 1)
      #average auc over all k
      cv_auc_sum <- cv_auc_sum + cv_auc
    }
    # substract current penalty from averaged auc
    cv_penalty <- penalty(a, w_candidates[j], balance, data)
    cv_auc_average <- cv_auc_sum/k - cv_penalty
    message("Current number of SNPs:", w_candidates[j],"  Average AUC:", cv_auc_average, " (Penalty:", cv_penalty, ")")

  }
  # get all aucs
  get_models <- sapply(X = 1:grid.res,
                       FUN = modelling)

  # select max auc
  ## welche spalte? annahme: 1. spalte w, 2. spalte auc
  idx_w <- max(get_models[,2])
  message("Best result for w = ",get_models[idx_w,1],"( AUC = ",get_models[idx_w,2],")")

  # save final model for prediction
  fin_model <- riskScorer(formula, data[, 1:w_candidates[idx_w]], y.name, feature.names, importance, weight, beta = TRUE)
  return(fin_model)

}