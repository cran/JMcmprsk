##' Joint modelling for longitutal and censored data with competing risks
##' @title ANOVA of longitudinal model
##' @param object  The JMcmprsk object returned by either jmo or jmc function.
 ##' @param coeff  Types of coefficients selected for anova. Note "alpha" is only avaiable to jmo type JMcmprsk object.
 ##' @param ... further arguments passed to or from other methods.
##' @return Return a Wald test statistic and the p value
##'   \tabular{ll}{
##'       \code{beta}    \tab  The Wald test for fixed effects for the longitutal part,i.e.  \eqn{\beta} in jmo or jmc output. \cr
##'       \code{gamma}    \tab  The Wald test for fixed effects for the survival part,i.e.  \eqn{\gamma} in jmo or jmc output. "gamma1" stands for test for competing risk 1 and "gamma2" stands for test for competing risk 2  \cr
##'       \code{alpha}    \tab  The Wald test for non-proportional odds covariates,i.e.  \eqn{\alpha} in jmo output. \cr
##'   }
##' @export
anova.JMcmprsk <-
  function (object,coeff=c("beta","gamma","alpha"), ...) {
    if (!inherits(object, "JMcmprsk"))
      stop("Use only with 'JMcmprsk' objects.\n")
    if (object$type == "jmo") {
      if (coeff=="beta") {
        betas <- object$betas
        nbeta=length(betas)
        betacov<- object$vcmatrix[1:nbeta,1:nbeta]
        W=t(betas)%*%(solve(betacov))%*%betas
       
        pval = pchisq(W,nbeta,lower.tail = FALSE)
        WaldTest <- data.frame(Chisq = W, df = nbeta, 
                          "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="beta" )
        

      }else if (coeff=="alpha") {
        nbeta=length(object$betas)
        #number of nonpropotional hazards 
        ncomp=nrow(object$alphamatrix)
        WaldTest<-NULL
        for(j in 1: ncomp){  
          alphas=as.vector(t(object$alphamatrix[j,]))
          nalpha=length(alphas)
          
          alphacov<- object$vcmatrix[(nbeta+1+(j-1)*nalpha):(nbeta+j*nalpha),(nbeta+1+(j-1)*nalpha):(nbeta+j*nalpha)]
          W=t(alphas)%*%(solve(alphacov))%*%alphas
          
          pval = pchisq(W,nalpha,lower.tail = FALSE)
          WaldTest <- rbind(WaldTest,data.frame(Chisq = W, df = nalpha, 
                                                "Pr(>|Chi|)" = pval, check.names = FALSE, row.names =paste0("alpha",j+1) ))
        }
        
        
      } else if (coeff=="gamma") {
        nbeta=length(object$betas)
        alphas=as.vector(t(object$alphamatrix))
        nalpha=length(alphas)
        #number of competing risks
        ncomp=nrow(object$gamma_matrix)
        WaldTest<-NULL
      for(j in 1: ncomp){  
        gammas=as.vector(t(object$gamma_matrix[j,]))
        ngamma=length(gammas)
        
        gammacov<- object$vcmatrix[(nbeta+nalpha+1+(j-1)*ngamma):(nbeta+nalpha+j*ngamma),(nbeta+nalpha+1+(j-1)*ngamma):(nbeta+nalpha+j*ngamma)]
        W=t(gammas)%*%(solve(gammacov))%*%gammas
        
        pval = pchisq(W,ngamma,lower.tail = FALSE)
        WaldTest <- rbind(WaldTest,data.frame(Chisq = W, df = ngamma, 
                               "Pr(>|Chi|)" = pval, check.names = FALSE, row.names =paste0("gamma",j) ))
      } 
  } else {
    print("Anova Error or Methods not developed yet!")
    WaldTest=NULL
   }      
      
      WaldTest
    } else if (object$type == "jmc")  {
      if (coeff=="beta") {
        #skip the first intercept
        betas <- object$betas
        nbeta=length(betas)-1
        betacov<- object$vcmatrix[2:nbeta,2:nbeta]
        W=t(betas[2:nbeta])%*%(solve(betacov))%*%betas[2:nbeta]
        
        pval = pchisq(W,nbeta,lower.tail = FALSE)
        WaldTest <- data.frame(Chisq = W, df = nbeta, 
                               "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="beta" )
        
        
      } else if (coeff=="gamma") {
        nbeta=length(object$betas)
   
        #number of competing risks
        ncomp=nrow(object$gamma_matrix)
        WaldTest<-NULL
        for(j in 1: ncomp){  
          gammas=as.vector(t(object$gamma_matrix[j,]))
          ngamma=length(gammas)
          
          gammacov<- object$vcmatrix[(nbeta+1+1+(j-1)*ngamma):(nbeta+1+j*ngamma),(nbeta+1+1+(j-1)*ngamma):(nbeta+1+j*ngamma)]
          W=t(gammas)%*%(solve(gammacov))%*%gammas
          
          pval = pchisq(W,ngamma,lower.tail = FALSE)
          WaldTest <- rbind(WaldTest,data.frame(Chisq = W, df = ngamma, 
                                                "Pr(>|Chi|)" = pval, check.names = FALSE, row.names =paste0("gamma",j) ))
        }
        
      }    else {
        print("Anova Error or Methods not developed yet!")
        WaldTest=NULL
      }     
      WaldTest
    }
}
