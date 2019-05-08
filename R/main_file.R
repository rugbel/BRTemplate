# display version number and date when the package is loaded
#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  BRTemplate\n',
    'Version:  ', desc$Version, '\n',
    'Date:     ', desc$Date, '\n',
    'Authors:  Ruggero Bellio (University of Udine)\n')
}


#' @keywords internal
splitdata <- function(X, y, ncores)
{
  M <- nrow(X)
  rest <- M %% ncores
  B.size <- floor(M / ncores)
  ind <-  if(ncores>1) sort(c(rep(1:(ncores-1), B.size), rep(ncores, B.size+rest)))
          else rep(1,B.size)
  listXy <- split(cbind(X, y), ind)
  return(listXy)
}




#' Generate the AD object corresponding to the model of interest
#'
#' @param dll Name of a  .cpp file  containing the TMB model template.
#' @param X Model matrix.
#' @param y Observed response.
#' @param parain Starting values (guesses) for parameter vector.
#' @export
get_AD <- function(cppfile, X, y, parain)
{
  lf <- nchar(cppfile)
  dll <- substr(cppfile,1,lf-4)
  compile(cppfile)
  dyn.load(dynlib(dll))
  parameters <- list(para = parain)
  data.list <- list(X=as.matrix(X), y=as.vector(y))
  obj <- MakeADFun(data=data.list, parameters=parameters, DLL=dll, silent=TRUE)
  return(obj)
}


#' Generates a single simulated data set merging all the MC data together
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @param y Observed response.
#' @param datagen Function that generates a single dataset. It must have two arguments,
#' the parameter value \code{para} and the model matrix \code{X}.
#' @param R Size of Monte Carlo simulation. Default is 500.
#' @export
gen_data <- function(para, X, y, datagen, R=500)
{
 Xout <- matrix(rep(t(X), R), ncol = ncol(X), byrow = TRUE)
 yout <- datagen(para, Xout)
 return(list(X=Xout, y=yout))
}


#' Computes the expected Fisher information
#'
#' Computes the expected Fisher information divided by the sample size.
#'
#' The function is just a portion of the \code{A_star} function.
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @param y Observed response.
#' @param ncores Number of cores used in estimation.
#' @param dll Name of the C++ template representing the model of interest.
#' @return The returned value is a vector of the same size of \code{para}, containing the
#' value of the adjustment.
#' @export
#' @author Ruggero Bellio
#' @references  Kosmidis, I. and Firth, D. (2010).  A generic algorithm for reducing bias
#' in parametric estimation.  \emph{ Electron. J. Statist.}, 4, 1097-1112.
#' @examples
#' # See the vignette file
exp_info1 <- function(para, X, y, ncores, dll)
{
  listXy <- splitdata(X, y, ncores)
  M <- nrow(X)
  p <- length(para)
  parameters <-  list(para = para)
  myf <- function(v) {
     matXy <- matrix(v, ncol=ncol(X)+1, byrow=FALSE)
     B.size <- nrow(matXy)
     Xv <- matXy[,1:ncol(X)]
     yv <- matXy[,ncol(X)+1]
     out1 <- matrix(0, B.size, p)
     for(i in 1:B.size){
       data.list1 <- list(X = matrix(Xv[i,], nrow=1), y = matrix(yv[i], nrow=1))
       uni1 <- MakeADFun(data=data.list1, parameters=parameters, DLL=dll, silent=TRUE)
       out1[i,] <- -uni1$gr(para)
     }
     list(out1=out1)
  }
  ogg <-  parallel::mclapply(listXy, myf, mc.cores = ncores)
  out1 <- matrix(0, M, p)
  startind <- 0
  for(i in 1:length(ogg)) {
    	B.size <- nrow(ogg[[i]]$out1)
        ind <- startind + 1:B.size
        startind <- startind + B.size
        out1[ind,] <- ogg[[i]]$out1
     }
  J <- cov(out1) * (M-1) / M
  return(J)
 }


#' @keywords internal
tr <- function(x) sum(diag(x))



#' Computes the BR adjustment
#'
#' Computes the BR adjustement based on the expected Fisher information.
#'
#' The function can be used for both the computation of the adjustment based on the empirical
#' approximation to the expected values and  the computation based on the Monte Carlo approximation.
#' In the first case, \code{X} and \code{y} would refer to observed data, while in the second
#' case they would refer to the simulated data. This is handled by the \code{gsolv} function.
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @param y Observed response.
#' @param ncores Number of cores used in estimation.
#' @param dll Name of the C++ template representing the model of interest.
#' @return The returned value is a vector of the same size of \code{para}, containing the
#' value of the adjustment.
#' @export
#' @author Ruggero Bellio
#' @references  Kosmidis, I. and Firth, D. (2010).  A generic algorithm for reducing bias
#' in parametric estimation.  \emph{ Electron. J. Statist.}, 4, 1097-1112.
#' @examples
#' # See the vignette file
A_star <- function(para, X, y, ncores, dll)
{
 listXy <- splitdata(X, y, ncores)
  M <- nrow(X)
  p <- length(para)
  parameters <-  list(para = para)
  a <- rep(0, p)
  myf <- function(v) {
     matXy <- matrix(v, ncol=ncol(X)+1, byrow=FALSE)
     B.size <- nrow(matXy)
     Xv <- matXy[,1:ncol(X)]
     yv <- matXy[,ncol(X)+1]
     E <- array(0, dim=c(p, p, p))
     out1 <- matrix(0, B.size, p)
     for(i in 1:B.size){
       data.list1 <- list(X = matrix(Xv[i,], nrow=1), y = matrix(yv[i], nrow=1))
       uni1 <- TMB::MakeADFun(data=data.list1, parameters=parameters, DLL=dll, silent=TRUE)
       out1[i,] <- -uni1$gr(para)
       out2 <- uni1$he(para)
       t1 <- tcrossprod(out1[i,])
       JS <- out2 - t1
       E <- E + outer(JS, out1[i,])
       }
     list(E=E, out1=out1)
    }
  ogg <-  parallel::mclapply(listXy, myf, mc.cores = ncores)
  startind <- 0
  E <- array(0, dim=c(p, p, p))
  out1 <- matrix(0, M, p)
  for(i in 1:length(ogg))
    {
        B.size <- nrow(ogg[[i]]$out1)
        ind <- startind + 1:B.size
        startind <- startind + B.size
        E <- E + ogg[[i]]$E / M
        out1[ind,] <- ogg[[i]]$out1
   }
  J <- cov(out1) * (M-1) / M
  F1 <- try(solve(J))
  if(is.numeric(F1)) for(j in 1:p) a[j] <- -0.5 * tr(F1 %*% E[j,,])
  return(a)
 }




#' Estimating function for BR estimation
#'
#' Estimating function for BR estimation following Kosmidis and Firth (2010).
#'
#' The function is written for usage with the \code{nleqslv} package. Both quasi Newton-Raphson
#' and quasi Fisher-scoring update are supported.
#'
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @param y Observed response.
#' @param ncores Number of cores to use. Recourse to the \code{parallel::detectCores} function if you're unsure
#' about them.
#' @param ADobj Object returned by \code{TMB::MakeADFun} based on the model specified by the
#' \code{dll} file. See the vignette file.
#' @param datagen Function that generates a single data set from the model of interest, based on the
#' the specified parameter value and the model matrix.
#' @param dll Name of the C++ template representing the model of interest.
#' @param R Number of simulated data sets for Monte Carlo computation of expected values. Not used if
#' \code{seed} is NULL.
#' @param seed Random seed for Monte Carlo computation. Default is NULL, which switches to the empirical
#' approximation.
#' @param observed Should the observed information used adopted for the direction of the update in the
#' estimation algorithm? If \code{FALSE}, the expected Fisher information is used instead. Default is
#' \code{TRUE}.
#' @return The estimating equation at \code{para}.
#' @export
#' @author Ruggero Bellio
#' @references  Kosmidis, I. and Firth, D. (2010).  A generic algorithm for reducing bias
#' in parametric estimation.  \emph{ Electron. J. Statist.}, 4, 1097-1112.
#' @examples
#' # See the vignette file
gsolv <- function(para, X, y, ncores, ADobj, dll, datagen, R, seed=NULL, observed=TRUE)
{
  if(!is.null(seed)){
         set.seed(seed)
         data_i <- gen_data(para, X, y, datagen, R=R)
         Xi <- data_i$X
         yi <- data_i$y
      } else {
              Xi <- X
              yi <- y
      }
  a <- A_star(para, Xi, yi, ncores, dll)
  score <- as.vector(-ADobj$gr(para))
  out <- score + a
  return(out)
}



#' Jacobian for the estimating function for BR estimation
#'
#' Jacobian for the estimating function for BR estimation.
#'
#' The function is written for usage with the \code{nleqslv} package. Both quasi Newton-Raphson
#' and quasi Fisher-scoring update are supported.
#'
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @param y Observed response.
#' @param ncores Number of cores to use. Recourse to the \code{parallel::detectCores} function if you're unsure
#' about them.
#' @param ADobj Object returned by \code{TMB::MakeADFun} based on the model specified by the
#' \code{dll} file. See the vignette file.
#' @param datagen Function that generates a single data set from the model of interest, based on the
#' the specified parameter value and the model matrix.
#' @param dll Name of the C++ template representing the model of interest.
#' @param R Number of simulated data sets for Monte Carlo computation of expected values. Not used if
#' \code{seed} is NULL.
#' @param seed Random seed for Monte Carlo computation. Default is NULL, which switches to the empirical
#' approximation.
#' @param observed Should the observed information used adopted for the direction of the update in the
#' estimation algorithm? If \code{FALSE}, the expected Fisher information is used instead. Default is
#' \code{TRUE}.
#' @return The estimating equation at \code{para}.
#' @export
#' @author Ruggero Bellio
#' @references  Kosmidis, I. and Firth, D. (2010).  A generic algorithm for reducing bias
#' in parametric estimation.  \emph{ Electron. J. Statist.}, 4, 1097-1112.
#' @examples
#' # See the vignette file
jsolv <- function(para, X, y, ncores, ADobj, dll, datagen, R, seed=NULL, observed = FALSE)
{
  if(!is.null(seed)){
         set.seed(seed)
         data_i <- gen_data(para, X, y, datagen, R=R)
         Xi <- data_i$X
         yi <- data_i$y
      } else {
              Xi <- X
              yi <- y
      }
  out <- if(!observed) exp_info1(para, Xi, yi, ncores, dll) * nrow(X)
         else ADobj$he(para)
  return(-out)
 }

