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
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
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
    tempj <- .Call("CCleanIndices", as.double(rows), as.double(nr), PACKAGE="bigmemory")
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
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
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
iptBM <- function(x)
{
  if(!is.big.matrix(x)){
    stop("Error: x must be of class 'big.matrix'")
  }
  
  # Call in-place-transpose
  .Call('CIPTMatrix', x@address)
  
  diff <- abs(ncol(x) - nrow(x))
  
  # resize big.matrix to transpose dimensions
  if(nrow(x) > ncol(x)){
    resizeBM(x, -diff, diff)
  }else{
    resizeBM(x, diff, -diff)
  }
}




#' @title Divide a "big.matrix" by an numeric value
#' @description This was written to provide division capabilities to big.matrix objects.  To avoid making another 
#' large matrix, this a wrapper of C++ around the big.matrix object to subtract each element from the provided value.  
#' It likely could use further optimization. It also includes column and row selection in addition to optionally filebacking.
#' @param x A \code{"big.matrix"}
#' @param value Numeric value (e.g. "1L")
#' @param cols Possible subset of columns for the transpose; could be numeric, named, or logical
#' @param rows Possible subset of rows for the transpose; could be numeric, named, or logical
#' @param y Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type preferably specified (e.g. "integer", "double", etc.)
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @return a \code{"big.matrix"} of type double
#' @export
divideIntBM <- function(x, value, cols=NULL, rows=NULL, 
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
    .Call("CDivisionIntBM", x@address, y@address, as.integer(value),
          getOption("bigmemory.typecast.warning"))
  else
    for (i in 1:length(cols)) y[,i] <- x[rows,cols[i]]
  
  return(y)
}

#' @title Divide a "big.matrix" by another "big.matrix"
#' @description This was written to provide division capabilities to big.matrix objects.  To avoid making another 
#' large matrix, this a wrapper of C++ around the big.matrix object to subtract each element from the provided value.  
#' It likely could use further optimization. It also includes column and row selection in addition to optionally filebacking.
#' @param x A \code{"big.matrix"}
#' @param y A \code{"big.matrix"}
#' @param x.cols Possible subset of columns for the transpose; could be numeric, named, or logical
#' @param x.rows Possible subset of rows for the transpose; could be numeric, named, or logical
#' @param y.cols Possible subset of columns for the transpose; could be numeric, named, or logical
#' @param y.rows Possible subset of rows for the transpose; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type preferably specified (e.g. "integer", "double", etc.)
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @return a \code{"big.matrix"} of type double
#' @export
divideMatrixBM <- function(x, y, x.cols=NULL, x.rows=NULL, 
                        y.cols=NULL, y.rows=NULL,
                        z=NULL, type=NULL, separated=NULL,
                        backingfile=NULL, backingpath=NULL,
                        descriptorfile=NULL, binarydescriptor=FALSE,
                        shared=TRUE)
{
  
  if(typeof(x) != "double" | typeof(y) != "double"){
    stop("Error: big.matrix types must be 'double'")
  }
  
  # subsetting matrices
  cols <- cleanupcols(x.cols, ncol(x), colnames(x))
  cols2 <- cleanupcols(y.cols, ncol(y), colnames(y))
  rows <- cleanuprows(x.rows, nrow(x), rownames(x))
  rows2 <- cleanuprows(y.rows, nrow(y), rownames(y))
  
  if(length(cols) != length(cols2)){
    stop("Error: number of columns not equal between matrices")
  }
  if(length(rows) != length(rows2)){
    stop("Error: number of rows not equal between matrices")
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
  
  # subsetting matrices
  if(length(cols) != ncol(x) | length(rows) != nrow(x)){
    xs <- deepcopy(x, cols=cols, rows=rows)
  }else{
    xs <- x
  }

  if(length(cols2) != ncol(y) | length(rows2) != nrow(y)){
    ys <- deepcopy(y, cols=cols2, rows=rows2)
  }else{
    ys <- y
  }

  
  if (is.null(z)) {
    z <- big.matrix(nrow=length(rows), ncol=length(cols), type=type, init=NULL,
                    dimnames=dimnames(x), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("CDivisionMatrixBM", xs@address, ys@address, z@address,
        getOption("bigmemory.typecast.warning"))
  
  return(z)
}


#' @title Multiply a "big.matrix" by another "big.matrix"
#' @description This was written to provide multiplication capabilities to big.matrix objects.  To avoid making another 
#' large matrix, this a wrapper of C++ around the big.matrix object to subtract each element from the provided value.  
#' It likely could use further optimization. It also includes column and row selection in addition to optionally filebacking.
#' @param x A \code{"big.matrix"}
#' @param y A \code{"big.matrix"}
#' @param x.cols Possible subset of columns for the transpose; could be numeric, named, or logical
#' @param x.rows Possible subset of rows for the transpose; could be numeric, named, or logical
#' @param y.cols Possible subset of columns for the transpose; could be numeric, named, or logical
#' @param y.rows Possible subset of rows for the transpose; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type preferably specified (e.g. "integer", "double", etc.)
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @return a \code{"big.matrix"} of type double
#' @export
multiplyMatrixBM <- function(x, y, x.cols=NULL, x.rows=NULL, 
                           y.cols=NULL, y.rows=NULL,
                           z=NULL, type=NULL, separated=NULL,
                           backingfile=NULL, backingpath=NULL,
                           descriptorfile=NULL, binarydescriptor=FALSE,
                           shared=TRUE)
{
  
  if(typeof(x) != "double" | typeof(y) != "double"){
    stop("Error: big.matrix types must be 'double'")
  }
  
  # subsetting matrices
  cols <- cleanupcols(x.cols, ncol(x), colnames(x))
  cols2 <- cleanupcols(y.cols, ncol(y), colnames(y))
  rows <- cleanuprows(x.rows, nrow(x), rownames(x))
  rows2 <- cleanuprows(y.rows, nrow(y), rownames(y))
  
  if(length(cols) != length(cols2)){
    stop("Error: number of columns not equal between matrices")
  }
  if(length(rows) != length(rows2)){
    stop("Error: number of rows not equal between matrices")
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
  
  # subsetting matrices
  if(length(cols) != ncol(x) | length(rows) != nrow(x)){
    xs <- deepcopy(x, cols=cols, rows=rows)
  }else{
    xs <- x
  }
  
  if(length(cols2) != ncol(y) | length(rows2) != nrow(y)){
    ys <- deepcopy(y, cols=cols2, rows=rows2)
  }else{
    ys <- y
  }
  
  
  if (is.null(z)) {
    z <- big.matrix(nrow=length(rows), ncol=length(cols), type=type, init=NULL,
                    dimnames=dimnames(x), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("CMultiplyMatrixBM", xs@address, ys@address, z@address,
        getOption("bigmemory.typecast.warning"))
  
  return(z)
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
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
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

#' @title cbind 'in-place' for class "big.matrix"
#' @description This is used to bind columns to an existing big.matrix object.  
#' WARNING!!! This modifies the original big.matrix object and therefore must be
#' either recreated or subset to regain the original structure.  Use only if you 
#' intend to only use the original matrix in the new form.
#' @param x A \code{"big.matrix"} 
#' @param y Columns to append, accepts \code{"big.matrix"} and \code{"matrix"} in addition to 
#' \code{"numeric"} and \code{"integer"} (e.g. "2.5", "1L")
#' @param cols.y Optionally specify columns to bind from 'y' if 'y' is a \code{"matrix"} or 
#' \code{"big.matrix"}.
#' @return The original \code{"big.matrix"} with the added columns
#' @example Examples/cbindBMIP.R
#' @export
cbindBMIP <- function(x, y, cols.y=NULL)
{
  if(class(x) != "big.matrix"){
    stop("x is not a 'big.matrix'")
  }
  
  if(is.big.matrix(y) | is.matrix(y)){
    if(nrow(x) != nrow(y)){
      stop("Error: Number of rows between x and y are not equal.")
    }
  }
  
  # set number of columns to add
  if(class(y) == "numeric" | class(y) == "integer"){
    newCols <- 1
  }else{
    if(!is.null(cols.y)){
      cols2 <- cleanupcols(cols.y, ncol(y), colnames(y))  
      newCols <- length(cols2)
    }else{
      newCols <- ncol(y)
      cols2 <- seq(newCols)
    }
  }
 
  # Add the needed number of columns
  resizeBM(x, 0, newCols)
  
  # Fill in added columns
  switch(class(y),
         integer = x[,ncol(x)] <- y,
         numeric = x[,ncol(x)] <- y,
         matrix = x[,(ncol(x)-newCols+1):ncol(x)] <- y[,c(cols2)],
         big.matrix = x[,(ncol(x)-newCols+1):ncol(x)] <- y[,c(cols2)],
         )
  
  # matrix modification message returned by 'resizeBM' so NULL return
}


#' @title Provide power (i.e. '^') for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to a specified power, with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param value A \code{"numeric"} value to take each element to the power
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type preferably specified (e.g. "integer", "double", etc.)
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been taken to the power
#' of \code{"value"}.  This function should only need to be used if the new matrix must be filebacked as the 
#' NULL call is applied in the \code{"^"} method.
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[bigmemoryExt]{powBMIP}} for in-place power
#' @export
powBM <- function(x, value,
                  cols.x=NULL,
                  z=NULL, type=NULL, separated=NULL,
                  backingfile=NULL, backingpath=NULL,
                  descriptorfile=NULL, binarydescriptor=FALSE,
                  shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  if (is.null(type)){
    type <- "double"
  } 
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }

  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }

  .Call("CpowBM", y@address, z@address, as.double(value), getOption("bigmemory.typecast.warning"))  
  
  return(z)
}


#' @title Provide in-place power (i.e. '^') for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to a specified power.
#' WARNING!!! This replaces the original matrix.
#' @param x A \code{"big.matrix"}
#' @param value A \code{"numeric"} value to take each element to the power
#' @details This function takes each element of a \code{"big.matrix"} to the specified power.  
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[base]{^}} for new big.matrix object \code{\link[bigmemoryExt]{powBM}} for more flexible new big.matrix
#' @export
powBMIP <- function(x, value)
{
  if(!is.big.matrix(x)){
    stop("Error: x is not of type 'big.matrix'")
  }
  
  if(!typeof(x) == "double"){
    stop("Error: the 'big.matrix' x must be of type 'double'")
  }
  
  if(!is.numeric(value) & !is.integer(value)){
    stop("Erorr: value must be of type 'integer' or 'numeric'")
  }
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  .Call("CpowBMIP", x@address, as.integer(value), getOption("bigmemory.typecast.warning"))  
  ret <- paste(c(paste("big.matrix", x@address, "was modified"), collapse=" "))
  return(ret)
}


#' @title Provide exponential function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} and return the exponential with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type default: "double"
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been replaced with its' exponent.  
#' This function should only need to be used if the new matrix must be filebacked as the 
#' NULL call is applied in the \code{\link[base]{exp}} method.
#' @return a \code{"big.matrix"}
#' @export
expBM <- function(x,
                  cols.x=NULL,
                  z=NULL, type="double", separated=NULL,
                  backingfile=NULL, backingpath=NULL,
                  descriptorfile=NULL, binarydescriptor=FALSE,
                  shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("CexpBM", y@address, z@address, as.double(seq(nrow(y))), as.double(cols1), getOption("bigmemory.typecast.warning"))  
  
  return(z)
}


#' @title Provide in-place exponential function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to a specified power.
#' WARNING!!! This replaces the original matrix.
#' @param x A \code{"big.matrix"} which must be of type "double" to allow decimals
#' @details This function takes each element of a \code{"big.matrix"} to its' respective exponential.  
#' @return a \code{"big.matrix"}
#' @export
expBMIP <- function(x)
{
  if(!is.big.matrix(x)){
    stop("Error: x is not of class 'big.matrix'")
  }
  
  if(typeof(x) != "double"){
    stop("Error: x must be of type 'double'")
  }

  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  .Call("CexpBMIP", x@address, getOption("bigmemory.typecast.warning"))  
  ret <- paste(c("big.matrix", x@address, "was modified"), collapse=" ")
  return(ret)
}

#' @title Provide common logarithm function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} and return the natural log with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type default: "double"
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been replaced with its' common log.  
#' This function should only need to be used if the new matrix must be filebacked as the NULL call is applied in the 
#' \code{\link[base]{log10}} method.  For other bases, use the \code{\link[bigmemoryExt]{logBaseBM}} function.
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[base]{log10}} for generic call, \code{\link[bigmemoryExt]{logBaseBM}} for more control of new matrix,
#' and \code{\link[bigmemoryExt]{logBMIP}} for in-place log transforms
#' @export
log10BM <- function(x,
                  cols.x=NULL,
                  z=NULL, type="double", separated=NULL,
                  backingfile=NULL, backingpath=NULL,
                  descriptorfile=NULL, binarydescriptor=FALSE,
                  shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("Clog10BM", y@address, z@address, getOption("bigmemory.typecast.warning"))  
  
  return(z)
}

#' @title Provide base logarithm function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} and return the natural log with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param base The base for the log calculation.  Default: \code{"exp(1)"}
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type default: "double"
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been replaced with its' natural log
#' of \code{"value"}.  This function should only need to be used if the new matrix must be filebacked as the 
#' NULL call is applied in the \code{\link[base]{log}} method.
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[base]{log10}} for generic common log10 call, \code{\link[bigmemoryExt]{log10BM}} for more control 
#' of new matrix, and \code{\link[bigmemoryExt]{logBMIP}} for in-place log transforms
#' @export
logBaseBM <- function(x, base,
                  cols.x=NULL,
                  z=NULL, type="double", separated=NULL,
                  backingfile=NULL, backingpath=NULL,
                  descriptorfile=NULL, binarydescriptor=FALSE,
                  shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }

  .Call("ClogBaseBM", y@address, z@address, base, getOption("bigmemory.typecast.warning"))  
  
  return(z)
}

#' @title Provide in-place natural log function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to a specified power.
#' WARNING!!! This replaces the original matrix.
#' @param x A \code{"big.matrix"} which must be of type "double" to allow decimals
#' @param base The base for the logarithm. Default = \code{"exp(1)"}
#' @details This function takes each element of a \code{"big.matrix"} to replaes with its' respective base logarithm.  
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[base]{log10}} for generic common log10 call, \code{\link[bigmemoryExt]{log10BM}} for more control 
#' of new matrix, and \code{\link[bigmemoryExt]{logBaseBM}} for different logartihm bases.
#' @export
logBMIP <- function(x, base=exp(1))
{
  if(!is.big.matrix(x)){
    stop("Error: x is not of class 'big.matrix'")
  }
  
  if(typeof(x) != "double"){
    stop("Error: x must be of type 'double'")
  }
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  .Call("ClogBMIP", x@address, base, getOption("bigmemory.typecast.warning"))  
  ret <- paste(c("big.matrix", x@address, "was modified"), collapse=" ")
  return(ret)
}


#' @title Provide hyperbolic tangent function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} and return the hyperbolic tangent with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type default: "double"
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been replaced with its' hyperbolic
#' tangent.  This function should only need to be used if the new matrix must be filebacked as the 
#' NULL call is applied in the \code{\link[base]{tanh}} method.
#' @return a \code{"big.matrix"}
#' @export
tanhBM <- function(x,
                  cols.x=NULL,
                  z=NULL, type="double", separated=NULL,
                  backingfile=NULL, backingpath=NULL,
                  descriptorfile=NULL, binarydescriptor=FALSE,
                  shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("CtanhBM", y@address, z@address, as.double(seq(nrow(y))), as.double(cols1), getOption("bigmemory.typecast.warning"))  
  
  return(z)
}


#' @title Provide in-place hyperbolic tangent function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to its' hyperbolic tangent.
#' WARNING!!! This replaces the original matrix.
#' @param x A \code{"big.matrix"} which must be of type "double" to allow decimals
#' @details This function takes each element of a \code{"big.matrix"} to its' hyperbolic tangent.  
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[bigmemoryExt]{tanhBM}} for separate big.matrix
#' @export
tanhBMIP <- function(x)
{
  if(!is.big.matrix(x)){
    stop("Error: x is not of class 'big.matrix'")
  }
  
  if(typeof(x) != "double"){
    stop("Error: x must be of type 'double'")
  }
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  .Call("CtanhBMIP", x@address, getOption("bigmemory.typecast.warning"))  
  ret <- paste(c("big.matrix", x@address, "was modified"), collapse=" ")
  return(ret)
}


#' @title Provide hyperbolic cosine function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} and return the hyperbolic cosine with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type default: "double"
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been replaced with its' hyperbolic
#' cosine.  This function should only need to be used if the new matrix must be filebacked as the 
#' NULL call is applied in the \code{\link[base]{cosh}} method.
#' @return a \code{"big.matrix"}
#' @export
coshBM <- function(x,
                   cols.x=NULL,
                   z=NULL, type="double", separated=NULL,
                   backingfile=NULL, backingpath=NULL,
                   descriptorfile=NULL, binarydescriptor=FALSE,
                   shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("CcoshBM", y@address, z@address, as.double(seq(nrow(y))), as.double(cols1), getOption("bigmemory.typecast.warning"))  
  
  return(z)
}


#' @title Provide in-place hyperbolic cosine function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to its' hyperbolic cosine.
#' WARNING!!! This replaces the original matrix.
#' @param x A \code{"big.matrix"} which must be of type "double" to allow decimals
#' @details This function takes each element of a \code{"big.matrix"} to its' hyperbolic cosine.  
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[bigmemoryExt]{tanhBM}} for separate big.matrix
#' @export
coshBMIP <- function(x)
{
  if(!is.big.matrix(x)){
    stop("Error: x is not of class 'big.matrix'")
  }
  
  if(typeof(x) != "double"){
    stop("Error: x must be of type 'double'")
  }
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  .Call("CcoshBMIP", x@address, getOption("bigmemory.typecast.warning"))  
  ret <- paste(c("big.matrix", x@address, "was modified"), collapse=" ")
  return(ret)
}


#' @title Provide hyperbolic sine function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} and return the hyperbolic sine with the 
#' option to have the new copy filebacked.
#' @param x A \code{"big.matrix"}
#' @param cols.x Possible subset of columns from the x object; could be numeric, named, or logical
#' @param z Optional destinitation object (matrix or big.matrix); if not specified, a big.matrix will be created
#' @param type default: "double"
#' @param separated use separated column organization of the data instead of column-major organization; 
#' use with caution if the number of columns is large.
#' @param backingfile the root name for the file(s) for the cache of x.
#' @param backingpath the path to the directory containing the file backing cache
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with 
#' \code{\link[bigmemory]{attach.big.matrix}}; if NULL, the backingfile is used as the root part of the descriptor file name. 
#' The descriptor file is placed in the same directory as the backing files.
#' @param binarydescriptor the flag to specify if the binary RDS format should be used for the backingfile description, 
#' for subsequent use with \code{\link[bigmemory]{attach.big.matrix}}; if NULL of FALSE, the dput() file format is used.
#' @param shared TRUE by default, and always TRUE if the big.matrix is file-backed. For a non-filebacked big.matrix, 
#' shared=FALSE uses non-shared memory, which can be more stable for large (say, >50% of RAM) objects. Shared memory 
#' allocation can sometimes fail in such cases due to exhausted shared-memory resources in the system.
#' @details This generates a new \code{"big.matrix"} object whereby each element has been replaced with its' hyperbolic
#' sine.  This function should only need to be used if the new matrix must be filebacked as the 
#' NULL call is applied in the \code{\link[base]{sinh}} method.
#' @return a \code{"big.matrix"}
#' @export
sinhBM <- function(x,
                   cols.x=NULL,
                   z=NULL, type="double", separated=NULL,
                   backingfile=NULL, backingpath=NULL,
                   descriptorfile=NULL, binarydescriptor=FALSE,
                   shared=TRUE)
{
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  if (is.null(type)) type <- typeof(x)
  if (is.big.matrix(x)) {
    if (is.null(separated)) separated <- is.separated(x)
  } else {
    separated <- FALSE
  }
  
  if(!is.null(cols.x)){
    cols1 <- cleanupcols(cols.x, ncol(x), colnames(x))
    y <- deepcopy(x, cols1)
  }else{
    cols1 <- seq(ncol(x))
    y <- x
  }
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(y), ncol=length(cols1), type=type, init=NULL,
                    dimnames=dimnames(y), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared)
  }
  
  .Call("CsinhBM", y@address, z@address, as.double(seq(nrow(y))), as.double(cols1), getOption("bigmemory.typecast.warning"))  
  
  return(z)
}


#' @title Provide in-place hyperbolic sine function for class "big.matrix"
#' @description This is used to take each element of a \code{"big.matrix"} to its' hyperbolic sine.
#' WARNING!!! This replaces the original matrix.
#' @param x A \code{"big.matrix"} which must be of type "double" to allow decimals
#' @details This function takes each element of a \code{"big.matrix"} to its' hyperbolic sine.  
#' @return a \code{"big.matrix"}
#' @seealso \code{\link[bigmemoryExt]{sinhBM}} for separate big.matrix
#' @export
sinhBMIP <- function(x)
{
  if(!is.big.matrix(x)){
    stop("Error: x is not of class 'big.matrix'")
  }
  
  if(typeof(x) != "double"){
    stop("Error: x must be of type 'double'")
  }
  
  if (nrow(x) > 2^31-1){
    stop(paste("Too many rows to copy at this point in time;",
               "this may be fixed in the future."))
  }
  
  .Call("CcoshBMIP", x@address, getOption("bigmemory.typecast.warning"))  
  ret <- paste(c("big.matrix", x@address, "was modified"), collapse=" ")
  return(ret)
}