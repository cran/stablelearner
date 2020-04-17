\name{similarity_values}
\alias{similarity_values}

\title{Extracting Similarity Values}

\description{
  Extract similarity values from object returned by \code{\link{stability}} for
  further illustration or analysis.
}

\usage{
  similarity_values(x, reverse = TRUE)
}

\arguments{
  \item{x}{an object of class \code{"stablelearner"} or 
    \code{"stablelearnerList"} from which the similarity values are extracted.}
  \item{reverse}{logical. If \code{reverse = TRUE} (default), the similarity 
    values are transformed (reversed) such that higher values indicate a higher 
    stability.}
}

\value{A numeric array of dimension 3 containing similarity values. The 
dimensions represent repetitions, results (fitted model objects) and similarity 
measures.}

\seealso{\code{\link{stability}}, \code{\link{summary.stablelearnerList}}}

\examples{
\donttest{
library("partykit")
res <- ctree(Species ~ ., data = iris)
stab <- stability(res)
similarity_values(stab)
}
}

\keyword{resampling}
\keyword{similarity}
