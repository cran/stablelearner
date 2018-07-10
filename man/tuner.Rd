\name{tuner}
\alias{tuner}

\title{Tuning Wrapper Function}

\description{
  Convenience function to train a method using different tuning parameters.
}

\usage{
  tuner(method, tunerange, \ldots)
}

\arguments{
  \item{method}{a character string. Name of the \R function to train the method.}
  \item{tunerange}{a list. A list that specifies the range of values to be used
	for each tuning parameter. Each element of the list should be a vector 
	that specifies the values to be tested for the tuning parameter. The 
	element must be named after the corresponding tuning parameter of the 
	method (see examples).}
  \item{\ldots}{additional information passed to \code{method} (such as
	\code{formula}, \code{data}, \code{subset}, etc.).}
}

\details{
  This function can be used to train any method using different values for its
  tuning parameter(s). The result can be passed directly to \code{\link{stability}}
  to compare the stability of results based on different values of the tuning 
  parameter.
}

\value{
  A list that contains all fitted model objects.

  Additional information about the range of values used for the tuning parameters
  is attached to the resulting object as an attribute.
}

\seealso{\code{\link{stability}}}

\examples{

\donttest{

library("partykit")

## tuning cforest using different values of its tuning parameter mtry
r <- tuner("cforest", tunerange = list(mtry = 1:4), formula = Species ~ ., data = iris)
stability(r, control = stab_control(seed = 1234))

## receive information about the range of tuning parameters
attr(r, "range")

}

}

\keyword{resampling, similarity}
