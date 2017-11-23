##' Print contents of JMcmprsk object.
##'
##'
##' @title Print JMcmprsk
##' @param x Object of class 'JMcmprsk'.
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{jmc}}
##' @author Hong Wang
##' @export
print.JMcmprsk <- function(x, ...) {
  if (!inherits(x, "JMcmprsk")) 
    stop("Not a legitimate \"JMcmprsk\" object")
  
  cat("Model Type:                            ", x$type, "\n\n")
  
  if (x$type == "jmc") {
    cat("                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
    cat("Longitudinal:                \n")
    cat(" Fixed effects:                 \n")
    
    for (i in 1:length(x$betas)) {
      #beta = paste0("beta", i)
      beta=names(x$betas)[i]
      uppsd = x$betas[i] + 1.96 * x$se_betas[i]
      lowersd = x$betas[i] - 1.96 * x$se_betas[i]
      zval = (x$betas[i]/x$se_betas[i])
      pval = 2 * pnorm(-abs(zval))
	  cat(" ",formatC(beta,width=14,flag="-"), sprintf("% 1.4f", x$betas[i]))	  
      cat("     ", sprintf("% 1.4f", x$se_betas[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    
    cat(" ",formatC("sigma^2",width=14,flag="-"), sprintf("% 1.4f", x$sigma2_val))
    cat("     ", sprintf("% 1.4f", x$se_sigma2_val))
    
    uppsd = x$sigma2_val + 1.96 * x$se_sigma2_val
    lowersd = x$sigma2_val - 1.96 * x$se_sigma2_val
    zval = (x$sigma2_val/x$se_sigma2_val)
    pval = 2 * pnorm(-abs(zval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("    ", pval)
    cat("\n")  
        
    
    # cat(' Estimate Std. Error 95%CI Pr(>|Z|) \n')
    cat("Survival:                \n")
    cat(" Fixed effects:                 \n")
    for (i in 1:dim(x$gamma_matrix)[1]) for (j in 1:dim(x$gamma_matrix)[2]) {
      gamma = paste0("gamma", i, j)
      gamma=colnames(x$gamma_matrix)[j]
      gamma=paste0(gamma,'_', i)
      gammaval = x$gamma_matrix[i, j]
      stdgammaval = x$se_gamma_matrix[i, j]
      uppsd = gammaval + 1.96 * stdgammaval
      lowersd = gammaval - 1.96 * stdgammaval
      zval = (gammaval/stdgammaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(gamma,width=14,flag="-"), sprintf("% 1.4f", gammaval))
      cat("     ", sprintf("% 1.4f", stdgammaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    cat("Random effects:                 \n")
	cat(" ",formatC("v2",width=14,flag="-"), sprintf("% 1.4f", x$v_estimate))
    cat("     ", sprintf("% 1.4f", x$se_v_estimate))
    
    uppsd = x$v_estimate + 1.96 * x$se_v_estimate
    lowersd = x$v_estimate - 1.96 * x$se_v_estimate
    zval = (x$v_estimate/x$se_v_estimate)
    pval = 2 * pnorm(-abs(zval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("     ", pval)
    cat("\n")
    
    # cat(' Random effects: \n')

    # restore a matrix from its uppertri
    sd_sigmamatrix = matrix(0, dim(x$sigma_matrix)[1], dim(x$sigma_matrix)[1])
    sd_sigmamatrix[lower.tri(sd_sigmamatrix, diag = TRUE)] <- x$se_sigma
    sd_sigmamatrix <- t(sd_sigmamatrix)
    
    # print sigmabii
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) # for (j in 1:dim(x$sigma_matrix)[2])
    {
      sigma = paste0("sigma_b", i, i)
      sigmaval = x$sigma_matrix[i, i]
      stdsigmaval = sd_sigmamatrix[i, i]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
	  cat(" ",formatC(sigma,width=14,flag="-"), sprintf("% 1.4f", sigmaval))
      cat("     ", sprintf("% 1.4f", stdsigmaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("    ", pval)
      cat("\n")
    }
    # print sigmau
    si = dim(x$sigma_matrix)[1]
    sigmaval = x$sigma_matrix[si, si]
    stdsigmaval = x$se_sigma[(si * (si + 1))/2]
    uppsd = sigmaval + 1.96 * stdsigmaval
    lowersd = sigmaval - 1.96 * stdsigmaval
    zval = (sigmaval/stdsigmaval)
    pval = 2 * pnorm(-abs(zval))
    cat(" ",formatC("sigma_u",width=14,flag="-"), sprintf("% 1.4f", sigmaval))
    cat("     ", sprintf("% 1.4f", stdsigmaval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("     ", pval)
    cat("\n")
    
    cat("Covariance:                 \n")
    
    
    # print sigmabii
    
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) for (j in (i + 
      1):(dim(x$sigma_matrix)[2])) {
      if (j < dim(x$sigma_matrix)[2]) {
        sigma = paste0("sigma_b", i, j)
      } else {
        sigma = paste0("sigma_b", i, "u")
      }
      
      sigmaval = x$sigma_matrix[i, j]
      stdsigmaval = sd_sigmamatrix[i, j]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(sigma,width=14,flag="-"), sprintf("% 1.4f", sigmaval))
      cat("     ", sprintf("% 1.4f", stdsigmaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
      
    }
    
  }
  if (x$type == "jmo") {    
    cat("                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
    cat("Longitudinal:                \n")
    cat(" Fixed effects:                 \n")
    
    for (i in 1:length(x$betas)) {
      #beta = paste0("beta", i)
      beta=names(x$betas)[i]
      uppsd = x$betas[i] + 1.96 * x$se_betas[i]
      lowersd = x$betas[i] - 1.96 * x$se_betas[i]
      zval = (x$betas[i]/x$se_betas[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(beta,width=14,flag="-"), sprintf("% 1.4f", x$betas[i]))	  
      cat("     ", sprintf("% 1.4f", x$se_betas[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
  
      for (i in 1:dim(x$alphamatrix)[1]) for (j in 1:dim(x$alphamatrix)[2]) {
      #alpha = paste0("alpha", i+1, j)
      alpha=colnames(x$alphamatrix)[j]
      alpha=paste0(alpha,'_', i+1)
      
      alphaval = x$alphamatrix[i, j]
      stdalphaval = x$se_alphas[(i-1)*dim(x$alphamatrix)[2]+j]
      uppsd = alphaval + 1.96 * stdalphaval
      lowersd = alphaval - 1.96 * stdalphaval
      zval = (alphaval/stdalphaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(alpha,width=14,flag="-"),sprintf("% 1.4f", alphaval))
      cat("     ", sprintf("% 1.4f", stdalphaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
  
  
  
  
  
	
	for (i in 1:length(x$thetas)) {
      theta = paste0("theta", i)
      uppsd = x$thetas[i] + 1.96 * x$se_thetas[i]
      lowersd = x$thetas[i] - 1.96 * x$se_thetas[i]
      zval = (x$thetas[i]/x$se_thetas[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(theta,width=14,flag="-"), sprintf("% 1.4f", x$thetas[i]))
      cat("     ", sprintf("% 1.4f", x$se_thetas[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    
    # cat(' Estimate Std. Error 95%CI Pr(>|Z|) \n')
    cat("Survival:                \n")
    cat(" Fixed effects:                 \n")
    for (i in 1:dim(x$gamma_matrix)[1]) for (j in 1:dim(x$gamma_matrix)[2]) {
      gamma = paste0("gamma", i, j)
      gamma=colnames(x$gamma_matrix)[j]
      gamma=paste0(gamma,'_', i)
      gammaval = x$gamma_matrix[i, j]
      stdgammaval = x$se_gamma_matrix[i, j]
      uppsd = gammaval + 1.96 * stdgammaval
      lowersd = gammaval - 1.96 * stdgammaval
      zval = (gammaval/stdgammaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(gamma,width=14,flag="-"), sprintf("% 1.4f", gammaval))
      cat("     ", sprintf("% 1.4f", stdgammaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    cat("Random effects:                 \n")
    cat(" ",formatC("v2",width=14,flag="-"), sprintf("% 1.4f", x$v_estimate))
    cat("      ", sprintf("% 1.4f", x$se_v_estimate))
    
    uppsd = x$v_estimate + 1.96 * x$se_v_estimate
    lowersd = x$v_estimate - 1.96 * x$se_v_estimate
    zval = (x$v_estimate/x$se_v_estimate)
    pval = 2 * pnorm(-abs(zval))
    cat("    ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("     ", pval)
    cat("\n")
    
    # cat(' Random effects: \n')
    
    # restore a matrix from its uppertri
    sd_sigmamatrix = matrix(0, dim(x$sigma_matrix)[1], dim(x$sigma_matrix)[1])
    sd_sigmamatrix[lower.tri(sd_sigmamatrix, diag = TRUE)] <- x$se_sigma
    sd_sigmamatrix <- t(sd_sigmamatrix)
    
    # print sigmabii
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) # for (j in 1:dim(x$sigma_matrix)[2])
    {
      sigma = paste0("sigma_b", i, i)
      sigmaval = x$sigma_matrix[i, i]
      stdsigmaval = sd_sigmamatrix[i, i]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(sigma,width=14,flag="-"), sprintf("% 1.4f", sigmaval))
      cat("      ", sprintf("% 1.4f", stdsigmaval))
      cat("    ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    # print sigmau
    si = dim(x$sigma_matrix)[1]
    sigmaval = x$sigma_matrix[si, si]
    stdsigmaval = x$se_sigma[(si * (si + 1))/2]
    uppsd = sigmaval + 1.96 * stdsigmaval
    lowersd = sigmaval - 1.96 * stdsigmaval
    zval = (sigmaval/stdsigmaval)
    pval = 2 * pnorm(-abs(zval))
    cat(" ",formatC("sigma_u",width=14,flag="-"), sprintf("% 1.4f", sigmaval))
    cat("      ", sprintf("% 1.4f", stdsigmaval))
    cat("    ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("     ", pval)
    cat("\n")
    
    cat("Covariance:                 \n")
    
    
    # print sigmabii
    
    for (i in 1:(dim(x$sigma_matrix)[1] - 1)) for (j in (i + 
                                                         1):(dim(x$sigma_matrix)[2])) {
      if (j < dim(x$sigma_matrix)[2]) {
        sigma = paste0("sigma_b", i, j)
      } else {
        sigma = paste0("sigma_b", i, "u")
      }
      
      sigmaval = x$sigma_matrix[i, j]
      stdsigmaval = sd_sigmamatrix[i, j]
      uppsd = sigmaval + 1.96 * stdsigmaval
      lowersd = sigmaval - 1.96 * stdsigmaval
      zval = (sigmaval/stdsigmaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ", formatC(sigma,width=14,flag="-"), sprintf("% 1.4f", sigmaval))
      cat("      ", sprintf("% 1.4f", stdsigmaval))
      cat("    ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
      
    }
    
    
  }
}

