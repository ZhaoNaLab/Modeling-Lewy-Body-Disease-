# subset columns of a matrix
# Mainly use for when we optionally subset the columns.
# In this case, the option subsample should be set to a number
subsample <- function(x, subsample=FALSE) {
  if (subsample) {
    keep = sample(1:ncol(x), size=subsample)
    if ('Seurat' %in% class(x)) {
      x <- subset(x, cells=colnames(x)[keep])
    } else {
      x <- x[,keep]
    }
  }
  return(x)
}

column_corr <- function(x, y=NULL, method=c('pearson', 'spearman')) {
  method = method[1]
  
  f <- function(a, b) {
    cA <- a - colMeans(a)
    cB <- b - colMeans(b)
    sA <- sqrt(colMeans(cA^2))
    sB <- sqrt(colMeans(cB^2))
    return(colMeans(cA * cB) / (sA * sB))
  }

  if (is.null(y)) { y = x } 
  if (method == 'spearman') {
    x = Rfast::colRanks(as.matrix(x))
    y = Rfast::colRanks(as.matrix(y))
  }
  return(f(x,y))
}

row_corr <- function(x, y=NULL, method=c('pearson', 'spearman')) {
  column_corr(x=x, y=y, method=method)
}

rowOrderStatistic <- function(x, n=1) {
  stopifnot(n > 0, n <= ncol(x))
  diag(x[,apply(x, MARGIN=1, function(z) { order(z)[n] })])
}

columnOrderStatistic <- function(x, n=1) {
  stopifnot(n > 0, n <= nrow(x))
  diag(x[apply(x, MARGIN=2, function(z) { order(z)[n] }),])
}
