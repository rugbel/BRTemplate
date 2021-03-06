% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_file.R
\name{gsolv}
\alias{gsolv}
\title{Estimating function for BR estimation}
\usage{
gsolv(para, X, y, ncores, ADobj, dll, datagen, R, seed = NULL,
  observed = TRUE)
}
\arguments{
\item{para}{Parameter vector.}

\item{X}{Model matrix.}

\item{y}{Observed response.}

\item{ncores}{Number of cores to use. Recourse to the \code{parallel::detectCores} function if you're unsure
about them.}

\item{ADobj}{Object returned by \code{TMB::MakeADFun} based on the model specified by the
\code{dll} file. See the vignette file.}

\item{dll}{Name of the C++ template representing the model of interest.}

\item{datagen}{Function that generates a single data set from the model of interest, based on the
the specified parameter value and the model matrix.}

\item{R}{Number of simulated data sets for Monte Carlo computation of expected values. Not used if
\code{seed} is NULL.}

\item{seed}{Random seed for Monte Carlo computation. Default is NULL, which switches to the empirical
approximation.}

\item{observed}{Should the observed information used adopted for the direction of the update in the
estimation algorithm? If \code{FALSE}, the expected Fisher information is used instead. Default is
\code{TRUE}.}
}
\value{
The estimating equation at \code{para}.
}
\description{
Estimating function for BR estimation following Kosmidis and Firth (2010).
}
\details{
The function is written for usage with the \code{nleqslv} package. Both quasi Newton-Raphson
and quasi Fisher-scoring update are supported.
}
\examples{
# See the vignette file
}
\references{
Kosmidis, I. and Firth, D. (2010).  A generic algorithm for reducing bias
in parametric estimation.  \emph{ Electron. J. Statist.}, 4, 1097-1112.
}
\author{
Ruggero Bellio
}
