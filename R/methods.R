
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