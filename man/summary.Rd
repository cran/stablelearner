\name{summary.stablelearnerList}
\alias{summary.stablelearner}
\alias{summary.stablelearnerList}
\alias{print.summary.stablelearnerList}

\title{Summarize Results from Stability Assessment}

\description{
  Summarizes and prints the results from stability assessments performed by 
  \code{\link{stability}}.
}

\usage{
\method{summary}{stablelearnerList}(object, ..., reverse = TRUE, 
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95), digits = 3, names = NULL)
}

\arguments{
  \item{object}{a object of class \code{"stablelearner"} or 
    \code{"stablelearnerList"} to be summarized.}
  \item{\dots}{Arguments passed from or to other functions (currently ignored).}
  \item{reverse}{logical. If \code{reverse = TRUE} (default), the similarity 
    values summarized in the \code{summary} output are transformed (reversed) 
    such that higher values indicate a higher stability.}
  \item{digits}{integer. The number of digits used to summarize the similarity 
    distribution in the \code{summary} output.}
  \item{probs}{a vector of probabilities used tosummarize the similarity 
    distribution in the \code{summary} output.}
  \item{names}{a vector of characters to specify a name for each result from 
    statistical learning in the \code{summary} output.}
}

\seealso{\code{\link{stability}}}

\examples{

\donttest{

library("partykit")

rval <- ctree(Species ~ ., data = iris)
stab <- stability(rval)

summary(stab)
summary(stab, reverse = FALSE)
summary(stab, probs = c(0.25, 0.5, 0.75))
summary(stab, names = "conditional inference tree")

}

}

\keyword{resampling}
\keyword{similarity}