
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

#' @export
exp.big.matrix <- function(e1){
  modArgs <- list(x=e1,
                  type="double")
  do.call("expBM", modArgs)
}

#' @export
log.big.matrix <- function(e1){
  modArgs <- list(x=e1,
                  type="double")
  do.call("logBM", modArgs)
}

