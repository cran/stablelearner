\name{dgp_twoclass}
\alias{dgp_twoclass}

\title{Data-Ggnerating Function for Two-Class Problem}

\description{
  Data-generating function to generate artificial data sets of a classification 
  problem with two response classes, denoted as \code{"A"} and \code{"B"}.
}

\usage{
  dgp_twoclass(n = 100, p = 4, noise = 16, rho = 0, 
    b0 = 0, b = rep(1, p), fx = identity)
}

\arguments{
  \item{n}{integer. Number of observations. The default is 100.}
  \item{p}{integer. Number of signal predictors. The default is 4.}
  \item{noise}{integer. Number of noise predictors. The default is 16.}
  \item{rho}{numeric value between -1 and 1 specifying the correlation 
    between the signal predictors. The correlation is given by \code{rho}^k, 
    where k is an integer value given by \code{\link{toeplitz}} 
    structure. The default is 0 (no correlation between predictors).}
  \item{b0}{numeric value. Baseline probability for class \code{"B"} on the logit 
    scale. The default is 0.}
  \item{b}{numeric value. Slope parameter for the predictors on the logit scale. 
    The default is 1 for all predictors.}
  \item{fx}{a function that is used to transform the predictors. The default
    is \code{\link{identity}} (equivalent to no transformation).}
}

\value{
  A \code{data.frame} including a column denoted as \code{class} that is
  a factor with two levels \code{"A"} and \code{"B"}. All other columns 
  represent the predictor variables (signal predictors followed by noise 
  predictors) and are named by \code{"x1"}, \code{"x2"}, etc..
}

\seealso{\code{\link{stability}}}

\examples{dgp_twoclass(n = 200, p = 6, noise = 4)}

\keyword{resampling}
\keyword{similarity}
