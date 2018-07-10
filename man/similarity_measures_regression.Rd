\name{similarity_measures_regression}

\alias{similarity_measures_regression}
\alias{edist}
\alias{msdist}
\alias{rmsdist}
\alias{madist}
\alias{qadist}
\alias{cprob}
\alias{ccc}
\alias{pcc}
\alias{cosine}
\alias{rbfkernel}
\alias{tanimoto}

\title{Similarity Measure Infrastructure for Stability Assessment with Numerical Responses}

\description{
  Functions that provide objects with functionality used by 
  \code{\link{stability}} to measure the similarity between numeric 
  predictions of two results in regression problems.
}

\usage{
  edist()
  msdist()
  rmsdist()
  madist()
  qadist(p = 0.95)
  
  cprob(kappa = 0.1)
  rbfkernel()
  tanimoto()
  cosine()
  
  ccc()
  pcc()
}

\arguments{
  \item{p}{A numeric value between 0 and 1 specifying the probability to which 
  the sample quantile of the absolute distance between the predictions is 
  computed.}
  \item{kappa}{A positive numeric value specifying the upper limit of the 
  absolute distance between the predictions to which the coverage probability 
  is computed.}
}

\details{

  The similarity measure functions provide objects that include functionality 
  used by \code{\link{stability}} to measure the similarity between numeric 
  predictions of two results in regression problems.

  The \code{edist} (euclidean distance), \code{msdist} (mean squared distance), 
  \code{rmsdist} (root mean squared distance), \code{madist} (mean absolute 
  distance) and \code{qadist} (quantile of absolute distance) functions 
  implement scale-variant distance measures that are unbounded.

  The \code{cprob} (coverage probability), \code{rbfkernel} (gaussian radial 
  basis function kernel), \code{tanimoto} (tanimoto coefficient) and 
  \code{cosine} (cosine similarity) functions implement scale-variant distance 
  measures that are bounded.

  The \code{ccc} (concordance correlation coefficient) and \code{pcc} (pearson
  correlation coefficient) functions implement scale-invariant distance 
  measures that are bounded between 0 and 1.

}

\seealso{\code{\link{stability}}}

\examples{

\donttest{

set.seed(0)

library("partykit")

airq <- subset(airquality, !is.na(Ozone))
m1 <- ctree(Ozone ~ ., data = airq[sample(1:nrow(airq), replace = TRUE),])
m2 <- ctree(Ozone ~ ., data = airq[sample(1:nrow(airq), replace = TRUE),])

p1 <- predict(m1)
p2 <- predict(m2)

## euclidean distance
m <- edist()
m$measure(p1, p2)

## mean squared distance
m <- msdist()
m$measure(p1, p2)

## root mean squared distance
m <- rmsdist()
m$measure(p1, p2)

## mean absolute istance
m <- madist()
m$measure(p1, p2)

## quantile of absolute distance
m <- qadist()
m$measure(p1, p2)

## coverage probability
m <- cprob()
m$measure(p1, p2)

## gaussian radial basis function kernel
m <- rbfkernel()
m$measure(p1, p2)

## tanimoto coefficient
m <- tanimoto()
m$measure(p1, p2)

## cosine similarity
m <- cosine()
m$measure(p1, p2)

## concordance correlation coefficient
m <- ccc()
m$measure(p1, p2)

## pearson correlation coefficient
m <- pcc()
m$measure(p1, p2)

}

}

\keyword{stability, similariy, measures}
