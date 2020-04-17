\name{similarity_measures_classification}

\alias{similarity_measures}
\alias{similarity_measures_classification}
\alias{clagree}
\alias{ckappa}
\alias{bdist}
\alias{tvdist}
\alias{hdist}
\alias{jsdiv}

\title{Similarity Measure Infrastructure for Stability Assessment with Ordinal Responses}

\description{
  Functions that provide objects with functionality used by 
  \code{\link{stability}} to measure the similarity between the predictions
  of two results in classification problems.
}

\usage{
  clagree()
  ckappa()

  bdist()
  tvdist()
  hdist()
  jsdiv(base = 2)
}

\arguments{
  \item{base}{A positive or complex number: the base with respect to which 
  logarithms are computed. Defaults to 2.}
}

\details{

  The similarity measure functions provide objects that include functionality 
  used by \code{\link{stability}} to measure the similarity between the 
  probability predictions of two results in classification problems.

  The \code{clagree} and \code{ckappa} functions provide an object that can be 
  used to assess the similarity based on the predicted classes of two results.
  The predicted classes are selected by the class with the highest probability.

  The \code{bdist} (Bhattacharayya distance), \code{tvdist} (Total variation 
  distance), \code{hdist} (Hellinger distance) and \code{jsdist} 
  (Jenson-Shannon divergence) functions provide an object that can be 
  used to assess the similarity based on the predicted class probabilities of 
  two results.
  
}

\seealso{\code{\link{stability}}}

\examples{

\donttest{

set.seed(0)

## build trees
library("partykit")
m1 <- ctree(Species ~ ., data = iris[sample(1:nrow(iris), replace = TRUE),])
m2 <- ctree(Species ~ ., data = iris[sample(1:nrow(iris), replace = TRUE),])

p1 <- predict(m1, type = "prob")
p2 <- predict(m2, type = "prob")

## class agreement
m <- clagree()
m$measure(p1, p2)

## cohen's kappa
m <- ckappa()
m$measure(p1, p2)

## bhattacharayya distance
m <- bdist()
m$measure(p1, p2)

## total variation distance
m <- tvdist()
m$measure(p1, p2)

## hellinger distance
m <- hdist()
m$measure(p1, p2)

## jenson-shannon divergence
m <- jsdiv()
m$measure(p1, p2)

## jenson-shannon divergence (base = exp(1))
m <- jsdiv(base = exp(1))
m$measure(p1, p2)

}

}

\keyword{stability}
\keyword{similariy}
\keyword{measures}
