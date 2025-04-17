#' Creates a list of row numbers for the V-fold cross validation.
#' @param N Sample size.
#' @param id Optional cluster id variable. If present, all observations in the same cluster will always be in the same split.
#' @param Y outcome.
#' @param cvControl Control parameters for the cross-validation step. See SuperLearner.CV.control for details.
#' 
#' @return A list of length V where each element in the list is a vector with the row numbers of the corresponding validation sample.
# created by Eric Polley on 2011-01-18.
CVFolds <- function (N, id, Y, cvControl) 
{
  if (!is.null(cvControl$validRows)) {
    return(cvControl$validRows)
  }
  stratifyCV <- cvControl$stratifyCV
  shuffle <- cvControl$shuffle
  V <- cvControl$V
  if (!stratifyCV) {
    if (shuffle) {
      if (is.null(id)) {
        validRows <- split(sample(1:N), rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(sample(1:n.id), rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
    else {
      if (is.null(id)) {
        validRows <- split(1:N, rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(1:n.id, rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
  }
  else {
    if (length(unique(Y)) != 2) {
      stop("stratifyCV only implemented for binary Y")
    }
    if (sum(Y) < V | sum(!Y) < V) {
      stop("number of (Y=1) or (Y=0) is less than the number of folds")
    }
    if (shuffle) {
      if (is.null(id)) {
        wiY0 <- which(Y == 0)
        wiY1 <- which(Y == 1)
        rowsY0 <- split(sample(wiY0), rep(1:V, length = length(wiY0)))
        rowsY1 <- split(sample(wiY1), rep(1:V, length = length(wiY1)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          validRows[[vv]] <- c(rowsY0[[vv]], rowsY1[[vv]])
        }
      }
      else {
        stop("stratified sampling with id not currently implemented")
      }
    }
    else {
      if (is.null(id)) {
        within.split <- suppressWarnings(tapply(1:N, 
                                                INDEX = Y, FUN = split, rep(1:V)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          validRows[[vv]] <- c(within.split[[1]][[vv]], 
                               within.split[[2]][[vv]])
        }
      }
      else {
        stop("stratified sampling with id not currently implemented")
      }
    }
  }
  return(validRows)
}

#' @param V Number of splits for the V-fold cross-validation step. The default is 10. In most cases, between 10 and 20 splits works well.
#' @param stratifyCV Should the data splits be stratified by a binary response? Attempts to maintain the same ratio in each training and validation sample.
#' @param shuffle	Should the rows of X be shuffled before creating the splits.
#' @param validRows Use this to pass pre-specified rows for the sample splits. The length of the list should be V and each entry in the list should contain a vector with the row numbers of the corresponding validation sample.
#' 
#' @return A list containing the control parameters.
# created by Eric Polley on 2011-01-18.
SuperLearner.CV.control<- function (V = 10L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL) 
{
  V <- as.integer(V)
  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("validRows must be a list of length V containing the row numbers for the corresponding validation set")
    }
    if (!identical(V, length(validRows))) {
      stop("V and length(validRows) must be identical")
    }
  }
  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}