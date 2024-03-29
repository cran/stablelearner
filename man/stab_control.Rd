\name{stab_control}
\alias{stab_control}

\title{Control for Supervised Stability Assessments}

\description{
  Various parameters that control aspects of the stability assessment performed
  via \code{\link{stability}}.
}

\usage{
  stab_control(B = 500, measure = list(tvdist, ccc), sampler = "bootstrap", 
    evaluate = "OOB", holdout = 0.25, seed = NULL, na.action = na.exclude,
    savepred = TRUE, silent = TRUE, ...)
}

\arguments{
  \item{B}{an integer value specifying the number of repetitions. The default
    is \code{B = 500}.}
  \item{measure}{a list of similarity measure (generating) functions. Those
    should either be functions of \code{p1} and \code{p2} (the predictions of 
    two results in a repetition, see Details below) or similarity measure 
    generator functions, see \code{\link{similarity_measures_classification}} 
    and \code{\link{similarity_measures_regression}}. The defaults are
    \code{\link{tvdist}} for the classification and \code{\link{ccc}} for the 
    regression case.}
  \item{sampler}{a resampling (generating) function. Either this should be a 
    function of \code{n} that returns a matrix or a sampler generator like 
    \code{\link{bootstrap}}. The default is \code{"bootstrap"}.}
  \item{evaluate}{a character specifying the evaluation strategy to be applied 
    (see Details below). The default is \code{"OOB"}.}
  \item{holdout}{a numeric value between zero and one that specifies the 
    proportion of observations hold out for evaluation over all repetitions,
    only if \code{evaluate = "OOB"}. The default is \code{0.25}.}
  \item{seed}{a single value, interpreted as an integer, see 
    \code{\link{set.seed}} and Details section below. The default is \code{NULL}, 
    which indicates that no seed has been set.}
  \item{na.action}{a function which indicates what should happen when the 
    predictions of the results contain \code{NAs}. The default function is
    \code{\link{na.exclude}}.}
  \item{savepred}{logical. Should the predictions from each iteration be 
    saved? If \code{savepred = FALSE}, the resulting object will be smaller.
    However, predictions are required to subsequently compute measures of 
    accuracy via \code{\link{accuracy}}.}
  \item{silent}{logical. If \code{TRUE}, error messages, generated by the 
    learner while training or predicting, are suppressed and \code{NA} is 
    returned for the corresponding iteration of the stability assessment 
    procedure.}
  \item{\dots}{arguments passed to \code{sampler} function.}
}

\details{

With the argument \code{measure} one or more measures can be defined that are
used to assess the stability of a result from supervised statistical learning
by \code{\link{stability}}. Predefined similarity measures for the regression
and the classification case are listed in \code{\link{similarity_measures_classification}} 
and \code{\link{similarity_measures_regression}}.

Users can define their own similarity functions \code{f(p1, p2)} that must 
return a single numeric value for the similarity between two results trained on 
resampled data sets. Such a function must take the arguments \code{p1} and \code{p2}. 
In the classification case, \code{p1} and \code{p2} are probability matrices of 
size m * K, where \code{m} is the number of predicted observations (size 
of the evaluation sample) and K is the number of classes. In the 
regression case, \code{p1} and \code{p2} are numeric vectors of length 
m.

A different way to implement new similarity functions for the current \R 
session is to define a similarity measure generator function, which is a
function without arguments that generates a list of five elements including the 
name of the similarity measure, the function to compute the similarity
between the predictions as described above, a vector of character values 
specifying the response types for which the similarity measure can be used, 
a list containing two numeric elements \code{lower} and \code{upper} that 
specify the range of values of the similarity measure and the function to 
invert (or reverse) the similarity values such that higher values indicate 
higher stability. The latter can be set to \code{NULL}, if higher similarity 
values already indicate higher stability. Those elements should be named
\code{name}, \code{measure}, \code{classes}, \code{range} and \code{reverse}.

The argument \code{evaluate} can be used to specify the evaluation strategy.
If set to \code{"ALL"}, all observations in the original data set are used for
evaluation. If set to \code{"OOB"}, only the pairwise out-of-bag observations
are used for evaluation within each repetition. If set to \code{"OOS"}, a 
fraction (defined by \code{holdout}) of the observations in the original data 
set are randomly sampled and used for evaluation, but not for training, over all 
repetitions.

The argument \code{seed} can be used to make similarity assessments comparable
when comparing the stability of different results that were trained on the same 
data set. By default, \code{seed} is set to \code{NULL} and the learning samples 
are sampled independently for each fitted model object passed to 
\code{\link{stability}}. If \code{seed} is set to a specific number, the seed
will be set for each fitted model object before the learning samples are 
generated using \code{"L'Ecuyer-CMRG"} (see \code{\link{set.seed}}) which 
guarantees identical learning samples for each stability assessment and, thus, 
comparability of the stability assessments between the results.}

\seealso{\code{\link{stability}}}

\examples{

\donttest{

library("partykit")
res <- ctree(Species ~ ., data = iris)

## less repetitions
stability(res, control = stab_control(B = 100))

}

\dontrun{

## change similarity measure
stability(res, control = stab_control(measure = list(bdist)))

## change evaluation strategy
stability(res, control = stab_control(evaluate = "ALL"))
stability(res, control = stab_control(evaluate = "OOS"))

## change resampling strategy to subsampling
stability(res, control = stab_control(sampler = subsampling))
stability(res, control = stab_control(sampler = subsampling, evaluate = "ALL"))
stability(res, control = stab_control(sampler = subsampling, evaluate = "OOS"))

## change resampling strategy to splithalf
stability(res, control = stab_control(sampler = splithalf, evaluate = "ALL"))
stability(res, control = stab_control(sampler = splithalf, evaluate = "OOS"))

}

}

\keyword{resampling}
\keyword{similarity}
