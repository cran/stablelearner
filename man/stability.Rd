\name{stability}
\alias{stability}
\alias{print.stablelearner}
\alias{print.stablelearnerList}

\title{Stability Assessment for Results from Supervised Statistical Learning}

\description{
  Stability assessment of results from supervised statistical learning (i.e., 
  recursive partitioning, support vector machines, neural networks, etc.). The 
  procedure involves the pairwise comparison of results generated from learning 
  samples randomly drawn from the original data set or directly from the 
  data-generating process (if available).
}

\usage{
  stability(x, \dots, data = NULL, control = stab_control(), weights = NULL, 
    applyfun = NULL, cores = NULL, names = NULL)
}

\arguments{
  \item{x}{fitted model object. Any model object can be used whose class is 
    registered in \code{\link{LearnerList}}. Users can add classes for the 
    current \R session to \code{\link{LearnerList}}, see \code{\link{addLearner}}.}
  \item{\ldots}{additional fitted model objects.}
  \item{data}{an optional \code{data.frame} or a data-generating function. By 
    default the learning data from \code{x} is used (if this can be inferred 
    from the \code{\link[stats]{getCall}} of \code{x}.}
  \item{control}{a list with control parameters, see  \code{\link{stab_control}}.}
  \item{weights}{an optional matrix of dimension n * B that can be used to 
    weight the observations from the original learning data when the models 
    are refitted. If \code{weights = true}, the weights are computed internally 
    according to the \code{sampler} defined in \code{control}. If 
    \code{weight = NULL} (default), no case-weights are used and the 
    \code{sampler} defined in \code{control} will be applied to the original 
    data set.}
  \item{applyfun}{a \code{\link[base]{lapply}}-like function. The default is to use
    \code{\link[base]{lapply}} unless \code{cores} is specified in which case
    \code{\link[parallel]{mclapply}} is used (for multicore computations on platforms
    that support these).}
  \item{cores}{integer. The number of cores to use in multicore computations
    using \code{\link[parallel]{mclapply}} (see above).}
  \item{names}{a vector of characters to specify a name for each fitted model 
    object. By default, the objects are named by their class.}
}

\details{
  FIXME
}

\value{
  For a single fitted model object, \code{stability} returns an object of 
  class \code{"stablelearner"} with the following components:
  
  \item{call}{the call from the model object \code{x},}
  \item{learner}{the information about the learner retrieved from \code{\link{LearnerList}},}
  \item{B}{the number of repetitions,}
  \item{sval}{a matrix containing the estimated similarity values for each 
  similarity measure specified in \code{control},}
  \item{sampstat}{a list containing information on the size of the learning 
  samples (\code{ls}), the size of the overlap between the learning samples 
  (\code{lo}), the size of the evaluation sample (\code{es}) and the size of
  the overlap between the evaluation and the learning samples (\code{eo})
  in each repetition.}
  \item{data}{a language object referring to the \code{data.frame} or the
    data-generating function used for assessing the stability,}
  \item{control}{a list with control parameters used for assessing the 
  stability,}
  
  For several fitted model objects, \code{stability} returns an object of 
  class \code{"stablelearnerList"} which is a list of objects of class 
  \code{"stablelearner"}.
}

\references{
Philipp M, Rusch T, Hornik K, Strobl C (2018). Measuring the
  \dQuote{Stability of Results from Supervised Statistical Learning}. Journal of
  Computational and Graphical Statistics. Forthcoming
}

\seealso{\code{\link{boxplot.stablelearnerList}}, \code{\link{summary.stablelearner}}}

\examples{

\donttest{
## assessing the stability of a single result
library("partykit")
r1 <- ctree(Species ~ ., data = iris)
stab <- stability(r1)
summary(stab)

## assessing the stability of several results
library("rpart")
r2 <- rpart(Species ~ ., data = iris)
stab <- stability(r1, r2, control = stab_control(seed = 0))
summary(stab, names = c("ctree", "rpart"))

## using case-weights instead of resampling
stability(r1, weights = TRUE)

## using self-defined case-weights
n <- nrow(iris)
B <- 500
w <- array(sample(c(0, 1), size = n*B*3, replace = TRUE), dim = c(n, B, 3))
stability(r1, weights = w)

## assessing stability for a given data-generating process
my_dgp <- function() dgp_twoclass(n = 100, p = 2, noise = 4, rho = 0.2)
res <- ctree(class ~ ., data = my_dgp())
stability(res, data = my_dgp)
}

}

\keyword{resampling, similarity}
