\name{boxplot.stablelearnerList}
\alias{boxplot.stablelearner}
\alias{boxplot.stablelearnerList}

\title{Illustrate Results from Stability Assessment}

\description{
  Illustrates the results from stability assessments performed by 
  \code{\link{stability}} using boxplots.
}

\usage{
  \method{boxplot}{stablelearnerList}(x, ..., main = NULL, xlab = NULL, ylab = NULL, reverse = TRUE)
  \method{boxplot}{stablelearner}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{"stablelearnerList"} to be illustrated.}
  \item{\dots}{Arguments passed to \code{\link{boxplot}}.}
  \item{main}{a character specifying the title. By default set to \code{NULL}.}
  \item{xlab}{a character specifying the title for the x axis. By default set 
    to \code{"Learner"}.}
  \item{ylab}{a character specifying the title for the y axis. By default set 
    to \code{NULL}.}
  \item{reverse}{logical. If \code{reverse = TRUE} (default), the similarity 
    values are transformed (reversed) such that higher values indicate a higher 
    stability.}
}

\seealso{\code{\link{stability}}, \code{\link{summary.stablelearnerList}}}

\examples{
\donttest{

library("partykit")
r1 <- ctree(Species ~ ., data = iris)

library("rpart")
r2 <- rpart(Species ~ ., data = iris)

stab <- stability(r1, r2, names = c("ctree", "rpart"))
boxplot(stab)

}
}

\keyword{resampling, similarity}
