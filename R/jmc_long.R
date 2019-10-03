#' An integrated function for reconstructing data and do the joint modelling 
#' @description Reconstruct data into a regular longitudinal format as a refined dataset and do joint modelling for this refined data with continuous outcome.  
#' @param long_data Data matrix for longitudinal in long form. The time variable should be labeled 'time'. 
#' @param surv_data Data matrix for competing risks data. Each subject has one row of observation (as opposed to the long_data).
#' First and second column should be the observed event time and censoring indicator, respectively. 
#' The coding for the censoring indicator is as follows: 0 - censored events, 1 - risk 1, 2 - risk 2. Two competing risks are assumed.
#' @param out Column name for outcome variable in long_data.
#' @param FE Vector of column names that correspond to the fixed effects in long_data. If missing, then all columns except for the outcome and ID columns will be considered.
#' @param RE Types/Vector of random effects in long_data. The available type are "intercept", "linear", "quadratic" (time-related random effect specification) or other covariates in the input dataset. If specify other covariates, then they to be numerical vectors.
#' @param ID Column name for subject ID number in long_data.
#' @param cate Vector of categorical variables in long_data. Default is NULL.
#' @param intcpt Specify either 0 or 1. Default is set as 1. 0 means no intercept in random effect.
#' @param quad.points Number of quadrature points used in the EM procedure. Default is 20. Must be an even number. Larger values means higher accuracy but more time-consuming.
#' @param max.iter Max iterations. Default is 10000.
#' @param quiet Logical. Print progress of function. Default is TRUE.
#' @param do.trace Logical. Print the parameter estimates during the iterations. Default is FALSE.
#' @return Object of class \code{JMcmprsk} with elements
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
##'@examples
##' \dontrun{
##' yfile=system.file("extdata", "rawfvc621_y.txt", package = "JMcmprsk")
##' cfile=system.file("extdata", "fvc621_c.txt", package = "JMcmprsk")
##' cread <- read.table(file = "rawfvc621_c.txt", header = T)
##' yread <- read.table(file = "fvc621_y.txt", header = T)
##' res <- jmc_long(long_data = yread, surv_data = cread, out = "FVC", cate = NULL,
##' FE = c("time", "FVC0", "FIB0", "CYC", "FVC0.CYC", "FIB0.CYC", "time.CYC"), 
##' RE = "time", ID = "rowId", intcpt = 1, quad.points = 8, max.iter = 10000, quiet = FALSE)
##' coef(res)
##' anova(res,coeff="beta")
##' anova(res,coeff="gamma")
##' #make up two categorical variables and add them into yread
##' require(tidyverse)
##' mfile=system.file("extdata", "fvc621_m.txt", package = "JMcmprsk")
##' mread <- read.table(file = "fvc621_m.txt", header = T)
##' rowId <- c(1:nrow(cread))
##' sex <- sample(c("Feamle", "Male"), nrow(mread), replace = T)
##' race <- sample(c("White", "Black", "Asian", "Hispanic"), nrow(mread), replace = T)
##' cate_var <- data.frame(rowId, sex, race)
##' yread <- left_join(yread, cate_var, by = "rowId")
##' # run jmc_long function again for yread file with two added categorical variables
##' res2 <- jmc_long(long_data = yread, surv_data = cread, out = "FVC", cate = c("sex", "race"), 
##' FE = c("time", "FVC0", "FIB0", "CYC", "FVC0.CYC", "FIB0.CYC", "time.CYC"), 
##' RE = "time", ID = "rowId", intcpt = 1, quad.points = 8, max.iter = 10000, quiet = TRUE)
##' }
##' @seealso \code{\link{jmc}}
##' @export

jmc_long <- function(long_data, surv_data, out, FE, RE, ID, cate = NULL,
                intcpt = 1, quad.points = 10, max.iter = 10000, quiet = TRUE, do.trace = FALSE) {
  
  if(!quiet) writeLines(">> Pre-processing data")
  # return the variable names for long_data
  long_names <- names(long_data)
  
  #Error Checking: if the typed variable couldn't be found in long_names, then the function will stop
  if(!out %in% long_names) 
    stop(paste0(out, " column is not included in long_data. Please correctly specify the outcome variable."))
  if(!is.numeric(long_data[, out])) 
    stop(paste0(out, "column should be numeric. Please correctly specify the outcome variable."))
  
  if(!is.null(FE)) {
    if(!prod(FE %in% long_names)) 
      stop(paste0(FE[!FE %in% long_names], " doesn't exist in long_data. Be sure fixed effects are correctly named."))
    for (j in FE)
    {
      if (!is.numeric(long_data[, j])) 
      {
        stop(paste0(j), " is not a numeric variable. Be sure to convert this variable to be numeric.")
        break
      }
    }
  }
  if(!RE %in% c(c("intercept", "linear", "quadratic"), long_names))
    stop(paste0(RE, " is an incorrect input for random effect setting. 
                Please choose intercept, linear, quadratic, or an existing variable(s) in long_data."))
  if(!ID %in% long_names) 
    stop(paste0(ID, " column is not included in long_data. Be sure ID variable is correctly named."))
  if(is.null(cate)) {
    writeLines(">> No categorical variables are considered in modelling...")
  } else if(!prod(cate %in% long_names)) {
    stop(paste0(cate[!cate %in% long_names], " column is not included in long_data. Be sure categorical variable is correctly named."))
  }
  if(!intcpt %in% c(0, 1))
    stop("intcpt argument can either be 0 or 1")
  if(max.iter <= 0)
    stop("max.iter must be a non-negative integer. Recommended is at least 1000.")
  if(quad.points %% 2 != 0)
    stop("quad.points must be an even number and not too small.")
  
  # create a matrix of fixed effects with a vector of 1's
  if(is.null(FE)) {
    fixed <- data.frame(1, long_data[, -which(colnames(long_data) %in% c(ID, out))])
  } else {
    fixed <- data.frame(1, long_data[, FE])
  }
  p     <- ncol(fixed)
  colnames(fixed)[1] <- "intercept"
  
  # create a matrix of random effects: generate new covariates based on the random effect specification
  if(RE == "intercept") {
    random <- data.frame(intcpt_RE = rep(1, nrow(long_data)))
  } else if(RE == "linear" & intcpt == 0) {
    if(!"time" %in% long_names)
      stop(paste0("To use RE = ", RE, "time must be named `time` in long_data"))
    random <- data.frame(time_RE = long_data[, "time"])
  } else if(RE == "linear" & intcpt == 1) {
    if(!"time" %in% long_names)
      stop(paste0("To use RE = ", RE, "time must be named `time` in long_data"))
    random <- data.frame(intcpt_RE = rep(1, nrow(long_data)), time_RE = long_data[, "time"])
  } else if(RE == "quadratic" & intcpt == 0) {
    if(!"time" %in% long_names)
      stop(paste0("To use RE = ", RE, "time must be named `time` in long_data"))
    random <- data.frame(time_RE = long_data[, "time"],
                         time2_RE = long_data[, "time"]^2)
  } else if(RE == "quadratic" & intcpt == 1) {
    if(!"time" %in% long_names)
      stop(paste0("To use RE = ", RE, "time must be named `time` in long_data"))
    random <- data.frame(intcpt_RE = rep(1, nrow(long_data)),
                         time_RE = long_data[, "time"],
                         time2_RE = long_data[, "time"]^2) 
  } else {
    random <- data.frame(long_data[, RE])
    colnames(random) <- paste0(RE, "_RE")
  }
  # the current version only allows for no more than 3 random effects
  if(ncol(random) > 3)
    stop("Current version can not have the number of effects exceeding 3.")
  
  #Format data
  ##Create dummy variables for the pre-determined categorical varibales
  ##Construct long_final with longitudinal outcome, fixed and random effects for jmc function
  ##Construct m_final for jmc function
  if(!is.null(cate)) {
    dummy_names <- c()
    for (j in cate){
      dummy <- matrix(0, nrow = nrow(long_data), ncol = (length(unique(long_data[, which(names(long_data) %in% j)]))-1))
      writeLines("The reference group for", j, "is", unique(long_data[, which(names(long_data) %in% j)])[1])
      for (i in 1:nrow(long_data))
      {
        if(as.numeric(long_data[, which(names(long_data) %in% j)])[i]>1){
          dummy[i, as.numeric(long_data[, which(names(long_data) %in% j)])[i]-1] <- 1 
        }
        else dummy[i, as.numeric(long_data[, which(names(long_data) %in% j)])[i]] <- 0 
      }
      colnames(dummy) <- sort(unique(long_data[, which(names(long_data) %in% j)]))[-1]
      dummy_names <- c(dummy_names, colnames(dummy))
      long_data <- cbind(long_data, dummy)
      long_data <- long_data[, -which(names(long_data) %in% j)]
    }
    outcome <- data.frame(long_data[, out])
    colnames(outcome) <- out
    dummy_variables <- data.frame(long_data[, dummy_names])
    p <- p+ncol(dummy_variables)
    long_final <- cbind(outcome, random, fixed, dummy_variables)
    m_final <- data.frame(n_count = table(long_data[, ID]))[, 2]
    m_final <- as.data.frame(m_final)
  } else {
    outcome <- data.frame(long_data[, out])
    colnames(outcome) <- out
    long_final <- cbind(outcome, random, fixed)
    m_final <- data.frame(n_count = table(long_data[, ID]))[, 2]
    m_final <- as.data.frame(m_final)
  }
  
  # check the survival times if they are all positive
  if(prod(surv_data[, 1] > 0) == 0) 
    stop("Survival time in surv_data (first column) must be positive valued.")
  
  if(!quiet) writeLines(">> Fitting the model")
  
  ##run jmc function
  res1 <- jmc(p = p, yfile = long_final, cfile = surv_data, 
                       mfile = m_final, point = quad.points, maxiterations = max.iter,
                       do.trace = do.trace, type_file = FALSE)
  
  ##return the joint modelling result
  return(res1)
}


