#' @useDynLib bigmemoryExt

cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc))
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}

cleanuprows <- function(rows=NULL, nr=NULL, rownames=NULL) {
  if (is.null(rows)) rows <- 1:nr
  else {
    if (!is.numeric(rows) & !is.character(rows) & !is.logical(rows))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(rows))
      if (is.null(rownames)) stop("row names do not exist.")
    else rows <- mmap(rows, rownames)
    if (is.logical(rows)) {
      if (length(rows) != nr)
        stop(paste("row vector length must match the number of",
                   "rows of the matrix."))
      rows <- which(rows)
    }
    tempj <- .Call("CCleanIndices", as.double(rows), as.double(nr))
    if (is.null(tempj[[1]])) stop("Illegal row index usage in extraction.\n")
    if (tempj[[1]]) rows <- tempj[[2]]
  }
  return(rows)
}


#' @title Produce transpose of a "big.matrix"
#' @description This is used to take the transpose of a big.matrix, with the new copy optionally filebacked
#' @param x A \code{"big.matrix"}
#' @param cols Possible subset of columns for the transpose; could be numeric, named, or logical
#' @param rows Possible subset of rows for the transpose; could be numeric, named, or logical
#' @param y Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type preferably specified (e.g. "integer", "double", etc.)
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath
#' @param descriptorfile
#' @param binarydescriptor
#' @param shared
#' @details This is need to make a tranpose of a \code{"big.matrix"} because R syntax won't access beyond the pointer.  
#' Converting between normal 'matrix' and 'big.matrix' also is expensive to use the normal \code{"t"} function.  This
#' also allows the user to optionally fileback the transposed object if it will be accessed several times.
#' @return a transposed \code{"big.matrix"}
#' @importMethodsFrom bigmemory typeof
#' @export
transposeBM <- function(x, cols=NULL, rows=NULL, 
                     y=NULL, type=NULL, separated=NULL,
                     backingfile=NULL, backingpath=NULL,
                     descriptorfile=NULL, binarydescriptor=FALSE,
                     shared=TRUE)
{
  cols <- cleanupcols(cols, ncol(x), colnames(x))
  rows <- cleanuprows(rows, nrow(x), rownames(x))
  if (nrow(x) > 2^31-1)
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  if (is.null(y)) {
    y <- big.matrix(nrow=length(cols), ncol=length(rows), type=type, init=NULL,
                    dimnames=dimnames(x), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  if (is.big.matrix(x) && is.big.matrix(y))
    .Call("CtransposeMatrix", x@address, y@address, as.double(rows), as.double(cols), 
          getOption("bigmemory.typecast.warning"))
  else
    for (i in 1:length(cols)) y[,i] <- x[rows,cols[i]]
  
  return(y)
}