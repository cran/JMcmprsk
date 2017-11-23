##' Joint modelling for longitutal and censored data with competing risks
##' @title Coefficients of longitudinal model
##' @param object  The JMcmprsk object returned by either jmo or jmc function.
##' @param ... further arguments passed to or from other methods.
##' @return Return estimates fixed effects with variable names
##' @export
coef.JMcmprsk <-
  function (object, ...) {
    if (!inherits(object, "JMcmprsk"))
      stop("Use only with 'JMcmprsk' objects.\n")
    if (object$type == "jmo") {
      betas <- object$betas
       out=betas
       out
    } else if (object$type == "jmc")  {
      betas <- object$betas
      out=betas
      out

    }
}
