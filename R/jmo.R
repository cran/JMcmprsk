##' Joint modeling of longitudinal ordinal data and competing risks
##' @title Joint Modelling for Ordinal outcomes 
##' @param p  The dimension of proportional odds covariates (not including intercept) in yfile.
##' @param s  The dimension of non-proportional odds covariates in yfile.
##' @param yfile Y matrix for longitudinal measurements in long format. For example, for a subject with n measurements, there are n rows for this subject. The # of rows in y matrix is the total number of measurements for all subjects. The columns in Y are ordered this way: the longitudinal outcome (column 1), then the covariates for random effects, and lastly, the covariates for fixed effects (no intercept).
##' @param cfile C matrix for competing risks failure time data. Each subject has one data entry, so the number of rows equals to the number of subjects. The survival / censoring time is included in the first column, and the failure type coded as 0 (censored events), 1 (risk 1), or 2 (risk 2) is given in the second column. Two competing risks are assumed. The covariates are included in the third column and on.
##' @param mfile M vector to indicate the number of longitudinal measurements per subject. The number of rows equals to the number of subjects.
##' @param point Quadrature points used in the EM procedure. Default is 20.
##' @param maxiterations Maximum values of iterations. Default is 100000.
##' @param do.trace Print detailed information of each iteration. Default is false, not to print the iteration details.
##' @param ... further arguments passed to or from other methods.
##' @return Object of class \code{JMcmprsk} with elements
##'   \tabular{ll}{
##'       \code{vcmatrix}    \tab  The variance-covariance matrix for all the parameters. The parameters are in the order: \eqn{\beta}, \eqn{\alpha}, \eqn{\theta}, \eqn{\gamma}, \eqn{\nu}, and \eqn{\Sigma}. The elements in \eqn{\Sigma} are output in the order along the main diagonal line, then the second main diagonal line, and so on. \cr
##'       \code{betas} \tab The point  estimates of \eqn{\beta}. \cr
##'       \code{se_betas} \tab The standard error estimate of \eqn{\beta}. \cr
##'       \code{alphamatrix} \tab The point  estimates of \eqn{\alpha}. \cr
##'       \code{se_alphas} \tab The standard error estimate of \eqn{\alpha}. \cr
##'       \code{theta} \tab The point  estimates of \eqn{\theta}. \cr 
##'       \code{se_theta} \tab The standard error estimate of \eqn{\theta}. \cr
##'       \code{gamma_matrix} \tab  The point  estimate of \eqn{\gamma}. \cr
##'       \code{se_gamma_matrix}   \tab  The standard error estimate of \eqn{\gamma}. \cr
##'       \code{v_estimate} \tab The point  estimate of \eqn{\nu}. \cr
##'       \code{se_v_estimate}    \tab The standard error estimate of \eqn{\nu}. \cr
##'       \code{sigma_matrix}     \tab The point estimate of \eqn{\Sigma} (only the upper triangle portion of the matrix is output).\cr
##'       \code{se_sigma}     \tab The standard error estimate of \eqn{\Sigma}.The standard errors are given in this order: main diagonal, the second main diagonal, and so on. \cr
##'       \code{loglike}     \tab Log Likelihood.\cr
##'   }
##' @examples
##' # A toy example on simulated data
##'  require(JMcmprsk)
##'  set.seed(123)
##'  yfile=system.file("extdata", "jmosimy.txt", package = "JMcmprsk")
##'  cfile=system.file("extdata", "jmosimc.txt", package = "JMcmprsk")
##'  mfile=system.file("extdata", "jmosimm.txt", package = "JMcmprsk")
##'  res3=jmo(p=3,s=1, yfile,cfile,mfile,point=10)
##'  coef(res3)
##'  anova(res3,coeff="beta")
##'  anova(res3,coeff="gamma")
##'  anova(res3,coeff="alpha")
##' #testing the function on real data with trace on
##'\dontrun{
##' require(JMcmprsk)
##' set.seed(123)
##' yfile=system.file("extdata", "ninds_nrank_y.txt", package = "JMcmprsk")
##' cfile=system.file("extdata", "ninds_nrank_c.txt", package = "JMcmprsk")
##' mfile=system.file("extdata", "ninds_nrank_m.txt", package = "JMcmprsk")
##' res=jmo(p=9,s=2, yfile,cfile,mfile,point=10,do.trace = TRUE)
##' coef(res)
##' anova(res,coeff="beta")
##' anova(res,coeff="gamma")
##' anova(res,coeff="alpha")
##' }
##' @references
##' \itemize{
##' \item Ning Li,Robert M. Elashoff,Gang Li and Jeffrey Saver. "Joint modeling of longitudinal ordinal data and competing risks survival times and analysis of the NINDS rt-PA stroke trial." Statistics in medicine 29.5 (2010): 546-557.
##' }
##' @seealso \code{\link{jmc}}
##' @export
jmo <- function (p,s,yfile,cfile,mfile,point=20,maxiterations=100000,do.trace=FALSE)
{
 
  # more error control here.
  # store header names for future useage
  ydata=read.table(yfile,header = T)
  ynames=colnames(ydata)
  yfile=tempfile(pattern = "", fileext = ".txt")
  writenh(ydata,yfile)
 
  cdata=read.table(cfile,header = T)
  cnames=colnames(cdata)
  cfile=tempfile(pattern = "", fileext = ".txt")
  writenh(cdata,cfile)

  ydim=dim(ydata)
  
  # number of observations in study is equals to the #of rows in Y matrix and delete header here
  n1=ydim[1];
  
  ydata[,1]=factor(ydata[,1])
  #generate data column names for further useage
  y_names=colnames(ydata)[1]
  fixed_col_names=paste(names(ydata[,(ncol(ydata)-p+1):ncol(ydata)]), collapse='+')
  initvalues=MASS::polr(as.formula(paste(y_names,"~",fixed_col_names)),data=ydata)
  betas=initvalues$coefficients
  thetas=initvalues$zeta
  
  # the levels of the first row of yfile
  K_num=length(unique(ydata[,1]))
  # dim of fixed effects plus dim of random effects should be 
  # total the column of y -the survival time column 1
  
  # dim of random effects
  p1a=ydim[2]-1-p-s;
  
  #if((p1<1)|(p1a<1)){
   # stop("Possibe wrong dimension of fixed effects in Y!")
  #}
  
  cdim=dim(cdata);
 
  # number of subjects in study is equals to the #of rows in C matrix
  k=cdim[1];
  # The dimension of fixed effects in C
  p2=cdim[2]-2;

  #The max number observations for a subject
  j_max=max(read.table(mfile));

 if (do.trace) { 
   trace=1;
 }else{
   trace=0;
 }
  
  
  myresult=jmo_main(k,n1, p,p2, p1a,s,K_num, j_max, point,betas,thetas, maxiterations, yfile,cfile,mfile,trace)
#the program is estimating -beta, -alpha, and -bi in equation (1) of the stats #in med paper.when reporting the results, the sign of beta, alpha (which is #beta2 in the program), and rho_bu (or sigma_bu, which is the off-diagonal  #elements of the sig matrix) should be flipped (i.e., negative values should be #positive, and vice versa)
  #beta and alpha
  myresult$betas=-myresult$betas
  myresult$alphamatrix=-myresult$alphamatrix
  
  #names
  names(myresult$betas)=ynames[(ydim[2]-p+1):ydim[2]]
  colnames(myresult$alphamatrix)=ynames[3:(s+3-1)]
  colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]
  
  
  #off-diagnoal elements
  
  for (i in 1:(dim(myresult$sigma_matrix)[1] - 1)) 
    for (j in (i +1):(dim(myresult$sigma_matrix)[2])) 
      {
      myresult$sigma_matrix[i,j]=-myresult$sigma_matrix[i,j]
    }
  
  
  myresult$type="jmo";
  
  class(myresult) <- "JMcmprsk"
  return (myresult)
}




