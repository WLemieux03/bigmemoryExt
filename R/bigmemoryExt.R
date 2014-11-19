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

#' @title 'Big.Matrix' In-Place-Transpose
#' @description This is used to apply and in-place transpose of a \code{"big.matrix"}.  It modifies the
#' original 'big.matrix' to its' corresponding tranpose.  WARNING!!! This changes your big.matrix object
#' and therefore the original structure is not accessible at the same time.  You can always convert back
#' to the original structure with this same function.
#' @param x A \code{"big.matrix"}
#' @param direction Specify which direction of transpose
#' @return Nothing, the \code{"big.matrix"} object has been modified.
#' @example Examples/iptExample.R
#' @export
iptBM <- function(x, direction)
{
  if(!is.big.matrix(x)){
    stop("Error: x must be of class 'big.matrix'")
  }
  if(!direction %in% c("r2c", "c2r")){
    stop("Error: 'direction' must by 'r2c' or 'c2r'.")
  }
  
  # Call in-place-transpose
  .Call('CIPTMatrix', x@address)
  
  # resize big.matrix to transpose dimensions
  switch(direction,
    c2r = resizeBM(x, -1, 1),
    r2c = resizeBM(x, 1, -1)
  )
#   if(direction == 'r2c'){
#     resizeBM(x, 1, -1)
#   }else if(direction == {
#     resizeBM(x, -1, 1) 
#   }
}

#' @title Subtract a "big.matrix" from an integer
#' @description This was written to provide subtraction capabilities to big.matrix objects.  The \code{"daxpy"}
#' function only works with equal size matrices.  To avoid making another large matrix, this a wrapper of C++ around the 
#' big.matrix object to subtract each element from the provided value.  It likely could use further optimization.
#' It also includes column and row selection in addition to optionally filebacking.
#' @param x A \code{"big.matrix"}
#' @param value A numeric value (e.g. "2.5", "1L")
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
#' @return a \code{"big.matrix"}
#' @export
subtractIntBM <- function(x, value, cols=NULL, rows=NULL, 
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
    y <- big.matrix(nrow=length(rows), ncol=length(cols), type=type, init=NULL,
                    dimnames=dimnames(x), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  if (is.big.matrix(x) && is.big.matrix(y))
    .Call("CsubtIntBM", x@address, y@address, as.double(rows), as.double(cols), as.integer(value),
          getOption("bigmemory.typecast.warning"))
  else
    for (i in 1:length(cols)) y[,i] <- x[rows,cols[i]]
  
  return(y)
}


#' @title cbind functionality for class "big.matrix"
#' @description This is used to bind columns to big.matrix objects, with the new copy optionally filebacked
#' @param x A \code{"big.matrix"}
#' @param y Columns to append, accepts \code{"big.matrix"} and \code{"matrix"} in addition to 
#' \code{"numeric"} and \code{"integer"} (e.g. "2.5", "1L")
#' @param binding Binding location, "right" or "left"
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param cols.y Possible subset of columns from the y object if a matrix type object otherwise ignored; 
#' could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type preferably specified (e.g. "integer", "double", etc.)
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath
#' @param descriptorfile
#' @param binarydescriptor
#' @param shared
#' @details This is simply a wrapper for \code{"big.matrix"} objects whereby a new \code{"big.matrix"} is
#' initialized and subsequently filled with the submitted arguments (optionally filtered).
#' @return a \code{"big.matrix"}
#' @export
cbindBM <- function(x, y, binding="right",
                    cols.x=NULL, cols.y=NULL,
                    z=NULL, type=NULL, separated=NULL,
                    backingfile=NULL, backingpath=NULL,
                    descriptorfile=NULL, binarydescriptor=FALSE,
                    shared=TRUE)
{
  
  if(class(y) == "big.matrix" | class(y) == "matrix"){
    if(nrow(x) != nrow(y)){
      stop(paste("The number of rows are not equal between matrices."))
    }
  }
  
  if (nrow(x) > 2^31-1)
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
  
  if(class(y) == "numeric" | class(y) == "integer"){
    cols2 <- 1
  }else{
    cols2 <- cleanupcols(cols.y, ncol(y), colnames(y))    
  }

  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x), ncol=length(cols1)+length(cols2), type=type, init=NULL,
                    dimnames=dimnames(x), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  if(binding=="right"){
    z[,1:length(cols1)] = x[,c(cols1)]
    
    if(is.big.matrix(y) | is.matrix(y)){
      z[,(length(cols1)+1):(length(cols1)+length(cols2))] = y[,c(cols2)]
    }
    else if(class(y) == "numeric" | class(y) == "integer"){
      z[,(length(cols1)+1)] = y
    }
  }else{
    if(binding == "left"){    
      if(is.big.matrix(y) | is.matrix(y)){
        z[,1:(length(cols2))] = y[,c(cols2)]
      }
      else if(class(y) == "numeric" | class(y) == "integer"){
        z[,1] = y
      }
      
      z[,(length(cols2)+1):(length(cols1)+length(cols2))] = x[,c(cols1)]
    }else{
      stop("binding variable not recognized, should be 'right' or 'left'")
    }
  }
  
  
  return(z)
}

# # Working here to create cbinding function in place
# #' @export
# cbindX <- function(x, y){
#   if(class(x) != "big.matrix"){
#     stop("x is not a 'big.matrix'")
#   }
#   if(class(y) != "matrix" & class(y) != "big.matrix"){
#     newCols <- 1
#   }else{
#     newCols <- ncol(y)
#   }
#  
#   resizeBM(x, nrow(x), newCols)
#   
#   x[,newCols] <- y
#   
# #   if(is.big.matrix(x) & is.big.matrix(y)){
# #     .Call("CBINDMatrix", x@address, y@address)
# #   }
# 
# }



