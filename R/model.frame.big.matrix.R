
#' @export
model.frame.big.matrix <- function(data,formula, drop.unused.levels = FALSE,na.action = na.fail,...){
  
  options(bigmemory.allow.dimnames=TRUE)
  
  possible_newdata <- !missing(data) && is.big.matrix(data) && 
    identical(deparse(substitute(data)), "newdata") && (nr <- nrow(data)) > 0
  if (!missing(formula) && nargs() == 1 && is.list(formula) && 
        !is.null(m <- formula$model)) 
    return(m)
  if (!missing(formula) && nargs() == 1 && is.list(formula) && 
        all(c("terms", "call") %in% names(formula))) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action"), names(fcall), 0)
    fcall <- fcall[c(1, m)]
    fcall[[1L]] <- quote(stats::model.frame)
    env <- environment(formula$terms)
    if (is.null(env)) 
      env <- parent.frame()
    return(eval(fcall, env))
  }
  
  if (missing(formula)){
    if (!missing(data) && inherits(data, "big.matrix") && 
          length(attr(data, "terms"))) 
      return(data)
    formula <- as.formula(data)
  }else{
    if (missing(data) && inherits(formula, "big.matrix")) {
      if (length(attr(formula, "terms"))) 
        return(formula)
      data <- formula
      formula <- as.formula(data)
    } 
  }
  formula <- as.formula(formula)
  
  if (missing(na.action)) {
    if (!is.null(naa <- attr(data, "na.action")) & mode(naa) != 
          "numeric") 
      na.action <- naa
    else if (!is.null(naa <- getOption("na.action"))) 
      na.action <- naa
  }
  if (missing(data)) {
    stop("you must provide a big.matrix object")
  }
  
  if (!inherits(formula, "terms")) 
    formula <- terms(formula)
  
  env <- environment(formula)
  rownames <- .row_names_info(data, 0L)
  vars <- attr(formula, "variables")
  predvars <- attr(formula, "predvars")
  if (is.null(predvars)) {
    predvars <- vars
  }
  
  varnames <- sapply(vars, function(x) paste(deparse(x, width.cutoff = 500), 
                                             collapse = " "))[-1L]
  
#   data <- as.big.matrix(data[, varnames])
  data <- deepcopy(data, cols=varnames)
  colnames(data) <- varnames
  return(data)
}