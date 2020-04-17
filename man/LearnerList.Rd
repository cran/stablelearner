\name{LearnerList}
\alias{LearnerList}

\title{List of Predefined Learners for Assessing Stability}

\description{
  The list contains details about several predefined learners that are required 
  to assess the stability of results from statistical learning.
}

\usage{LearnerList}

\details{
Currently implemented learners are:
\describe{
  \item{ctree}{conditional inference trees using \code{\link{ctree}} from \pkg{partykit}.}
  \item{rpart}{recursive partitioning using \code{rpart} from \pkg{rpart}.}
  \item{J48}{recursive partitioning using \code{J48} from \pkg{RWeka}.}
  \item{C5.0}{recursive partitioning using \code{C5.0} from \pkg{C50}.}
  \item{tree}{recursive partitioning using \code{tree} from \pkg{tree}.}
  \item{lda}{linear discriminant analysis using \code{lda} from \pkg{MASS}.}
  \item{lm}{linear models using \code{lm} from \pkg{stats}.}
  \item{glm}{generalized linear models using \code{glm} from \pkg{stats}.}
  }
Users can add new learners to \code{LearnerList} for the current \R session,
see \code{\link{addLearner}}.
}

\seealso{\code{\link{addLearner}}}

\keyword{resampling}
\keyword{similarity}
