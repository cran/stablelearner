\name{addLearner}
\alias{addLearner}

\title{Add Learners to \code{\link{LearnerList}}}

\description{
  The function can be used to add new learner to \code{\link{LearnerList}} in the 
  current \R session.
}

\usage{addLearner(x)}

\arguments{
  \item{x}{a list containing all required information to define a new learner 
    (see Details below).}
}

\details{
  The function can be used to add new learners to \code{\link{LearnerList}} in
  the current \R session. The function expects a list of four elements
  including the name of the learners object class, the name of the package 
  where the class and the fitting method is implemented, the name of the method 
  and a prediction function that predicts class probabilities (in the 
  classification case) or numeric values (in the regression case) and takes the 
  arguments \code{x} (the fitted model object), \code{newdata} a 
  \code{data.frame} containing the predictions of the observations in the
  evaluation sample and \code{yclass} a character string specifying the type
  of the response variable (\code{"numeric"}, \code{"factor"}, etc.). The 
  elements in the list should be named \code{class}, \code{package}, 
  \code{method} and \code{predfun}.
}

\seealso{\code{\link{LearnerList}}}

\examples{

newlearner <- list(
  class   = "svm",
  package = "e1071",
  method  = "Support Vector Machine",
  predict = function(x, newdata, yclass = NULL) {
    if(match(yclass, c("ordered", "factor"))) {
      attr(predict(x, newdata = newdata, probability = TRUE), "probabilities")
    } else {
      predict(x, newdata = newdata)
    }
  })

addLearner(newlearner)

}

\keyword{stability}
\keyword{similarity}
\keyword{resampling}
