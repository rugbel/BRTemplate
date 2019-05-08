#' Generates a single simulated data set from a logistic regression model
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @export
genLogit <- function(para, X)
{
 prob <- plogis(X %*% para)
 yout <- rbinom(nrow(X), size=1, prob=prob)
 return(yout)
}



#' Generates a single simulated data set from a probit regression model
#'
#' @param para Parameter vector.
#' @param X Model matrix.
#' @export
genProbit <- function(para, X)
{
  prob <- pnorm(X %*% para)
  yout <- rbinom(nrow(X), size=1, prob=prob)
  return(yout)
}
