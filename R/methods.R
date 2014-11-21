
##################
### S4 Methods ###
##################

#' @export
setMethod("Arith", c(e1="big.matrix", e2="numeric"),
          function(e1,e2)
            {
            op = .Generic[[1]]
            switch(op,
                   `^` = powBM(e1, e2),
                    stop("Undefined operation")
                   )
          }
          )

##################
### S3 Methods ###
##################

#####################################
### Exponential and Log Functions ###
#####################################

#' @export
exp.big.matrix <- function(e1){
  modArgs <- list(x=e1,
                  type="double")
  do.call("expBM", modArgs)
}

#' @export
log.big.matrix <- function(e1, base = exp(1)){
    modArgs <- list(x=e1,
                    base = base,
                    type="double")
    do.call("logBaseBM", modArgs)
}

#' @export
log10.big.matrix <- function(e1){
  modArgs <- list(x=e1,
                  type="double")
  do.call("log10BM", modArgs)
}

############################
### Hyperbolic Functions ###
############################

#' @export
tanh.big.matrix <- function(e1){
  tanhBM(e1)
}

#' @export
cosh.big.matrix <- function(e1){
  coshBM(e1)
}

#' @export
sinh.big.matrix <- function(e1){
  sinhBM(e1)
}

