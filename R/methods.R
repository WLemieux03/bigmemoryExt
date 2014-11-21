
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

#' @export
setMethod("Math", c(x="big.matrix"),
          function(x)
          {
            op = .Generic[[1]]
            switch(op,
                   `sqrt` = powBM(x, 0.5),
                   `exp` = expBM(x, type="double"),
                   `log10` = log10BM(x, type="double"),
                   #`log` = logBaseBM(x, type="double", base=base),
                   `tanh` = tanhBM(x),
                   `cosh` = coshBM(x),
                   `sinh` = sinhBM(x),
                   stop("Undefined operation")
            )
          }
)

#' @export
setMethod("log", signature(x="big.matrix"),
          function(x, base=exp(1))
            {
            logBaseBM(x, type="double", base=base)
          }
          )


