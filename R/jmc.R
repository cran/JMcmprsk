##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes 
##' @param p  The dimension of fixed effects (not including intercept) in yfile.
##' @param yfile Y matrix for longitudinal measurements in long format. For example, for a subject with n measurements, there should be n rows for this subject. The # of rows in y matrix is the total number of measurements for all subjects in the study. The columns in Y should start with the longitudinal outcome (column 1), the covariates for the random effects, and then the covariates for the fixed effects.
##' @param cfile C matrix for competing risks failure time data. Each subject has one data entry, so the number of rows equals to the number of subjects. The survival / censoring time is included in the first column, and the failure type coded as 0 (censored events), 1 (risk 1), or 2 (risk 2) is given in the second column. Two competing risks are assumed. The covariates are included in the third column and on.
##' @param mfile M vector to indicate the number of longitudinal measurements per subject. The number of rows equals to the number of subjects.
##' @param point Quadrature points used in the EM procedure.Default is 20.
##' @param maxiterations Maximum values of iterations. Default is 100000.
##' @param do.trace Print detailed information of each iteration. Default is false, i.e., not to print the iteration details.
##' @param type_file Types of inputs. Default is true, i.e.  data files with headers. If set to "F", inputs are changed to data matrixes or data.frames (with headers)
##' @return Object of class \code{JMcmprsk} with elements
##'   \tabular{ll}{
##'       \code{vcmatrix}    \tab  The variance-covariance matrix for all the parameters. The parameters are in the order: \eqn{\beta}, \eqn{\sigma^2}, \eqn{\gamma}, \eqn{\nu}, and \eqn{\Sigma}. The elements in \eqn{\Sigma} are output in the order along the main diagonal line, then the second main diagonal line, and so on. \cr
##'       \code{betas} \tab The point  estimates of \eqn{\beta}. \cr
##'       \code{se_betas} \tab The standard error estimate of \eqn{\beta}. \cr
##'       \code{gamma_matrix} \tab  The point  estimate of \eqn{\gamma}. \cr
##'       \code{se_gamma_matrix}   \tab  The standard error estimate of \eqn{\gamma}. \cr
##'       \code{v_estimate} \tab The point  estimate of \eqn{\nu}. \cr
##'       \code{se_v_estimate}    \tab The standard error estimate of \eqn{\nu}. \cr
##'       \code{sigma2_val}     \tab  The point estimate of \eqn{\sigma^2}.\cr
##'       \code{se_sigma2_val}     \tab  The standard error estimate of \eqn{\sigma^2}.\cr
##'       \code{sigma_matrix}     \tab The point estimate of \eqn{\Sigma} (only the upper triangle portion of the matrix is output).\cr
##'       \code{se_sigma}     \tab The standard error estimate of \eqn{\Sigma}.The standard errors are given in this order: main diagonal, the second main diagonal, and so on. \cr
##'       \code{loglike}     \tab Log Likelihood.\cr
##'   }
##'   
##' @examples
##' # A toy example on simulated data
##' require(JMcmprsk)
##' set.seed(123)
##' yfile=system.file("extdata", "jmcsimy.txt", package = "JMcmprsk")
##' cfile=system.file("extdata", "jmcsimc.txt", package = "JMcmprsk")
##' mfile=system.file("extdata", "jmcsimm.txt", package = "JMcmprsk")
##' res2=jmc(p=4,yfile,cfile,mfile,point=6)
##' coef(res2)
##' anova(res2,coeff="beta")
##' anova(res2,coeff="gamma")
##' #testing the function on real data with trace on
##' \donttest{
##' require(JMcmprsk)
##' set.seed(123)
##' yfile=system.file("extdata", "fvc621_y.txt", package = "JMcmprsk")
##' cfile=system.file("extdata", "fvc621_c.txt", package = "JMcmprsk")
##' mfile=system.file("extdata", "fvc621_m.txt", package = "JMcmprsk")
##' res1=jmc(p=7,yfile,cfile,mfile,do.trace = TRUE)
##' #if the input are not files but matrixes or data.frames,i.e. type_file=F
##'  ydata=read.table(yfile,header = T)
##'  cdata=read.table(cfile,header = T)
##'  mdata=read.table(mfile)
##'  res1=jmc(p=7,ydata,cdata,mdata, do.trace = TRUE,type_file = F)
##' coef(res1) 
##' anova(res1,coeff="beta") 
##' anova(res1,coeff="gamma")   
##' } 
##' @references
##' \itemize{
##' \item Elashoff, Robert M., Gang Li, and Ning Li. "A joint model for longitudinal measurements and survival data in the presence of multiple failure types." Biometrics 64.3 (2008): 762-771.
##' }
##' @seealso \code{\link{jmo}}
##' @export
jmc<- function (p,yfile,cfile,mfile,point=20,maxiterations=100000,do.trace=FALSE,type_file=TRUE)
{
  if (do.trace) { 
    trace=1;
  }else{
    trace=0;
  }
  
  
  #Gaussian-Hermite quadrature nodes and weights
  #The dimension of xs/ws is half of the point value since they are symmetric
  
  gq_vals <- statmod::gauss.quad(n = point, kind = "hermite")
  
  xs <- gq_vals$nodes[(point / 2 + 1) : point]
  
  ws <- gq_vals$weights[(point / 2 + 1) : point]
  
  if (type_file){
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
  # number of observations in study is equals to the #of rows in Y matrix
  n1=ydim[1];
  # dim of fixed effects plus dim of random effects should be 
  # total the column of y -the survival time column 1
  p1a=ydim[2]-1-p;
  
  if((p<1)|(p1a<1)){
    stop("Possibe wrong dimension of fixed effects in Y!")
  }
  
  cdim=dim(cdata);
 
  # number of subjects in study is equals to the #of rows in C matrix
  k=cdim[1];
  #
  p2=cdim[2]-2;

  maxl=max(read.table(mfile));
  myresult=jmc_main(k,n1, p,p2, maxl,p1a, maxiterations, point,xs,ws, yfile,cfile,mfile,trace)
  
  }else{
    ynames=colnames(yfile)
    yfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(yfile,yfilenew)
    
    cnames=colnames(cfile)
    cfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(cfile,cfilenew)
    
    mfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(mfile,mfilenew)
    
    ydim=dim(yfile)
    # number of observations in study is equals to the #of rows in Y matrix
    n1=ydim[1];
    # dim of fixed effects plus dim of random effects should be 
    # total the column of y -the survival time column 1
    p1a=ydim[2]-1-p;
    
    if((p<1)|(p1a<1)){
      stop("Possibe wrong dimension of fixed effects in Y!")
    }
    
    cdim=dim(cfile);
    
    # number of subjects in study is equals to the #of rows in C matrix
    k=cdim[1];
    #
    p2=cdim[2]-2;
    
    maxl=max(mfile);
  
    
  myresult=jmc_main(k,n1, p,p2, maxl,p1a, maxiterations, point,xs,ws, yfilenew,cfilenew,mfilenew,trace)  
  }
  
  

 


  myresult$type="jmc";

  #names
  names(myresult$betas)=ynames[(ydim[2]-p+1):ydim[2]]

  colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]

  myresult$k=k
  class(myresult) <- "JMcmprsk"
  
  return (myresult)
}




