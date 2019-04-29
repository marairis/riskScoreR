
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
  if(charmatch("x", feature.names) != "nomatch") {
    colnames(data) <- c(paste("x", (1:(ncol(data)-1)), sep = ""), y.name)
    names(importance) <-  paste("x", (1:(ncol(data)-1)), sep = "")
  } else if (charmatch("snp", feature.names) != "nomatch") {
    colnames(data) <- c(paste("snp", (1:(ncol(data)-1)), sep = ""), y.name)
    names(importance) <-  paste("snp", (1:(ncol(data)-1)), sep = "")
  } 
  
  # initialize penalty function
  penalty <- function(a, w, balance, data) {
    step <- ceiling(log(((w/a)+2)/2, 1.5))
    current_penalty <- step*balance
    return(current_penalty)
  }

  # (ii) initilize optimization problem (min)
  resp_matrix <- matrix(0, nrow = nrow(data), ncol = length(importance))
  
  dev_modelling <- function(w) {
    # true/false matrix, size depending on #snps = w
    tf_matrix <- matrix(1, nrow = w, ncol = w, byrow = TRUE)
    tf_matrix[lower.tri(tf_matrix)] <- 0
    # importance with size depending on w
    dev_importance <- sort(importance, decreasing = TRUE)[w]
    # data with size depending on w
    dev_data <- dplyr::select(data, names(dev_importance), y.name)
  
    # multiply true/false matrix with importance
    temp <- tf_matrix*dev_importance
  
    # response matrix
    resp_model <- riskScorer(data = dev_data,
                                   y.name = "state",
                                   feature.names = names(dev_importance),
                                   importance = dev_importance,
                                   weight = TRUE,
                                   beta = TRUE)
    
    # save response in matrix // einfach mit trainingsdaten predicten? 
    resp_values <- predict.riskScorer(resp_model,
                                      dev_data[1:(length(dev_importance))],
                                      "response")
    
    resp_matrix[,w] <- unlist(resp_values)
  }
  
  # matrix with all possible responses(w) for all person (n)
  dev_models <- sapply(X = 1:length(importance), 
                       FUN = dev_modelling)
  
  # get devianz for each model (#models = w) 
  # nur wert fuer family: gaussian/laplace (egal welcher type predictet wurde) // werte zwischen 0-1 
  # family: verteilung der zielvariable (siehe plot?)
  mod_deviance <- matrix(data = sapply(X = 1:(ncol(dev_models)-1),
                                       FUN = function (s) dismo::calc.deviance(obs = as.numeric(dev_models[,ncol(dev_models)]), 
                                                                               pred = as.numeric(dev_models[,s]), 
                                                                               family = "gaussian")),
                         ncol = ncol(dev_models)-1,
                         nrow = 1,
                         byrow = TRUE)
  
  # get penalty [auf 0-1 beschraenkt]
  max_pen <- penalty(a = a, w = length(importance), balance = balance, data = data)
  mod_penalty <- matrix(data = sapply(X = 1: ncol(mod_deviance),
                                      FUN = function (s) penalty(w = s, a = a, balance = balance, data = data)/max_pen),
                        ncol = ncol(mod_deviance),
                        nrow = 1,
                        byrow = TRUE)
  
  # solve optimization problem: min (devianz + penalty)
  best_model <- which.min(mod_deviance + mod_penalty)
  
  # get final model for prediction
  fin_importance <- sort(importance, decreasing = TRUE)[best_model]
  fin_feature.names <- names(fin_importance) 
  fin_y.name <- y.name
  fin_data <- dplyr::select(data, names(fin_importance), fin_y.name)
  
  fin_model <- riskScorer(data = fin_data, 
                          y.name = fin_y.name, 
                          feature.names = names(fin_importance),
                          importance = fin_importance,
                          weight = weight, 
                          beta = TRUE)
  return(fin_model)
}