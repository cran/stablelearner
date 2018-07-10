\name{accuracy}
\alias{accuracy}

\title{Prediction Accuracy from Stability Assessment Results}

\description{
  Function to compute the prediction accuracy from an object
  of class \code{"stablelearner"} or \code{"stablelearnerList"} as a parallel 
  to the similarity values estimated by \code{\link{stability}} in each 
  iteration of the stability assessment procedure.
}

\usage{
  accuracy(x, measure = "kappa", na.action = na.exclude, 
	   applyfun = NULL, cores = NULL) 
}

\arguments{
  \item{x}{an object of class \code{"stablelearner"} or 
	\code{"stablelearnerList"}.}
  \item{measure}{a character string (or a vector of character strings). 
	Name(s) of the measure(s) used to compute accuracy. Currently implemented
    	measures are \code{"diag"} = percentage of observations on the main diagonal 
	of a confusion matrix, \code{"kappa"} = \code{"diag"} corrected for agreement 
        by chance (default), \code{"rand"} = Rand index, and \code{"crand"} = 
	Rand index corrected for agreemend by chance 
	(see also \code{\link{classAgreement}}).}
  \item{na.action}{a function which indicates what should happen to the predictions 
	of the results containing \code{NAs}. The default function is 
	\code{\link{na.exclude}}.}
  \item{applyfun}{a \code{\link[base]{lapply}}-like function. The default is to use
    \code{\link[base]{lapply}} unless \code{cores} is specified in which case
    \code{\link[parallel]{mclapply}} is used (for multicore computations on platforms
    that support these).}
  \item{cores}{integer. The number of cores to use in multicore computations
    using \code{\link[parallel]{mclapply}} (see above).}
}

\details{
  This function can be used to compute prediction accuracy after the stability was
  estimated using \code{\link{stability}}.
}

\value{
  A matrix of size \code{2*B} times length(\code{measure}) containing prediction
  accuracy values of the learners trained during the stability assessment procedure.
}

\seealso{\code{\link{stability}}}

\examples{

\donttest{
library("partykit")
res <- ctree(Species ~ ., data = iris)
stab <- stability(res)
accuracy(stab)
}

}

\keyword{resampling, similarity}
