\name{bootstrap}
\alias{bootstrap}
\alias{samplesplitting}
\alias{jackknife}
\alias{subsampling}
\alias{splithalf}

\title{Sampler Infrastructure for Stability Assessment}

\description{
  Sampler objects that provide objects with functionality used by 
  \code{\link{stabletree}} to generate resampled datasets.
}

\usage{
  bootstrap(B = 500, v = 1)
  subsampling(B = 500, v = 0.632)
  samplesplitting(k = 5)
  jackknife(d = 1, maxrep = 5000)
  splithalf(B = 500)
}

\arguments{
  \item{B}{An integer value specifying the number of resampled datasets.}
  \item{k}{An integer value specifying the number of folds in sample-splitting.}
  \item{d}{An integer value specifying the number of observations left out in jackknife.}
  \item{maxrep}{An integer value specifying the maximum number of resampled datasets allowed, when using jackknife.}
  \item{v}{A numeric value between 0 and 1 specifying the fraction of observations in each subsample.}
}

\details{

  The sampler functions provide objects that include functionality to generate 
  resampled datasets used by \code{\link{stabletree}}.

  The \code{bootstrap} function provides an object that can be used to generate
  \code{B} bootstrap samples by sampling from \code{n} observations with 
  replacement.

  The \code{subsampling} function provides an object that can be used to 
  generate \code{B} subsamples by sampling from \code{floor(v*n)} 
  observations without replacement.

  The \code{samplesplitting} function provides an object that can be used to 
  generate \code{k}-folds from \code{n} observations.

  The \code{jackknife} function provides an object that can be used to generate
  all datasets necessary to perform leave-\code{k}-out jackknife sampling from 
  \code{n} observations. The number of datasets is limited by \code{maxrep} to 
  prevent unintended CPU or memory overload by accidently choosing too large 
  values for \code{k}.

  The \code{splithalf} function provides an object that can be used to 
  generate \code{B} subsamples by sampling from \code{floor(0.5*n)} 
  observations without replacement. When used to implement the "splithalf" 
  resampling strategy for measuring the stability of a result via the 
  \code{\link{stability}} function, the matrix containing the complement
  learning samples is generated automatically by \code{\link{stability}}.

}

\seealso{\code{\link{stabletree}}, \code{\link{stability}}}

\examples{

set.seed(0)

## bootstrap sampler
s <- bootstrap(3)
s$sampler(10)

## subsampling
s <- subsampling(3, v = 0.6)
s$sampler(10)

## 5-fold sample-splitting
s <- samplesplitting(5)
s$sampler(10)

## jackknife
s <- jackknife(d = 1)
s$sampler(10)

## splithaf
s <- splithalf(3)
s$sampler(10)

}

\keyword{regression}
