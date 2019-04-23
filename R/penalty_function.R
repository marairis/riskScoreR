
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
  checkmate::assertNumeric(importance, any.missing = FALSE)
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
  # make sure that feature.names and importance/training/testing column names match
  ## mit "x1" nicht ideal, allgemeine abfrage fuer x/snp?
  if("x1" %in% feature.names) {
    colnames(data) <- c(paste("x", (1:(ncol(data)-1)), sep = ""), y.name)
    names(importance) <-  paste("x", (1:(ncol(data)-1)), sep = "")
  } else if ("snp1" %in% feature.names) {
    colnames(data) <- c(paste("snp", (1:(ncol(data)-1)), sep = ""), y.name)
    names(importance) <-  paste("snp", (1:(ncol(data)-1)), sep = "")
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
  # calculate possible values for w
  w_candidates <- sapply(list(1:grid.res), function(s) ((ncol(data)-1)/grid.res)*s)
  cv_auc_penalized <- double(grid.res)
  
  # describe modelling function: does a cross-validation and returns performance for current w
  modelling <- function(j) {
    cv_auc_values <- double(k)
    for (i in 0:(k-1)) {
      # select train and test data
      cv_test <- data[perm_idx[((n_data/k)*i+1):((n_data/k)*(i+1))], ]
      cv_train <- data[-perm_idx[((n_data/k)*i+1):((n_data/k)*(i+1))], ]

      # set up model
      # create parameters adjusted to w
      #cv_formula <- as.formula(paste("state ~", paste(paste0("x", 1:w_candidates[j]), collapse = "+")))
      cv_y.name <- y.name
      ### in var.ranked: x statt snp, da spaltenname x in importance?
      cv_feature.names <- names(variables_ranked[1:w_candidates[j]])
      cv_importance <- variables_ranked[1:w_candidates[j]]
      cv_data <- dplyr::select(cv_train, cv_feature.names, y.name)
      
      cv_model <- riskScorer(data = cv_data,
                             y.name = cv_y.name,
                             feature.names = cv_feature.names,
                             importance = cv_importance,
                             weight = weight,
                             beta = TRUE)

      # validate on test data
      cv_newdata <- dplyr::select(cv_test, cv_feature.names, y.name)
      cv_prediction <- predict.riskScorer(risk.scorer = cv_model,
                                          newdata = cv_newdata,
                                          type = "prob")

      # get auc
      ## welche spalte uebergeben?
      cv_auc_values[(i+1)] <- mlr::measureAUC(probabilities = cv_prediction$`1`,
                                              truth = cv_test[, y.name],
                                              negative = 1,
                                              positive = 0)
    }
    
    # substract current penalty from averaged auc
    ## wenn auc < 0.5: 1-auc, dh modell umdrehen, wie?
    cv_penalty <- penalty(a = a,
                          w = w_candidates[j], 
                          balance = balance, 
                          data = data)
    cv_auc_penalized[j] <- mean(cv_auc_values) - cv_penalty
    if (cv_auc_penalized[j] < 0) {
      cv_auc_penalized[j] <- 0
    } 
    message("Current number of SNPs:", w_candidates[j],"  Averaged AUC:", cv_auc_penalized[j], " (Penalty:", cv_penalty, ")")
    return(cv_auc_penalized[j])
  }
  # get all aucs as an vector (vlt besser dataframe? mit anzahl snps, auc, penalty?)
  get_models <- sapply(X = 1:grid.res,
                       FUN = modelling)
  
  # noch machen: auc > 1, ausgabe finales modell -> auswahl aus get_models

  # select max auc
  idx_w <- which.max(get_models)
  message("Best result for w = " ,w_candidates[idx_w], " (AUC = ", get_models[idx_w], ")")

  # get final model for prediction
  fin_feature.names <- names(variables_ranked[1:w_candidates[idx_w]]) 
  fin_y.name <- y.name
  fin_importance <- variables_ranked[1:w_candidates[idx_w]]
  fin_data <- dplyr::select(data, fin_feature.names, fin_y.name)
  
  fin_model <- riskScorer(data = fin_data, 
                          y.name = fin_y.name, 
                          feature.names = fin_feature.names,
                          importance = fin_importance,
                          weight = weight, 
                          beta = TRUE)
  return(fin_model)
}