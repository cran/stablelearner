\name{stabletree}
\alias{stabletree}
\alias{print.stabletree}
\alias{summary.stabletree}
\alias{print.summary.stabletree}

\title{Stability Assessment for Tree Learners}

\description{
  Stability assessment of variable and cutpoint selection in
  tree learners (i.e., recursive partitioning). By refitting
  trees to resampled versions of the learning data, the stability
  of the trees structure is assessed and can be summarized and
  visualized.
}

\usage{
  stabletree(x, data = NULL, sampler = bootstrap, weights = NULL,
    applyfun = NULL, cores = NULL, savetrees = FALSE, \dots)
}

\arguments{
  \item{x}{fitted model object. Any tree-based model object that can be coerced
    by \code{\link[partykit]{as.party}} can be used provided that suitable methods
    are provided.}
  \item{data}{an optional \code{data.frame}. By default the learning data from \code{x}
    is used (if this can be inferred from the \code{\link[stats]{getCall}} of \code{x}.}
  \item{sampler}{a resampling (generating) function. Either this should be a function
    of \code{n} that returns a matrix or a sampler generator like \code{\link{bootstrap}}.}
  \item{weights}{an optional matrix of dimension n * B that can be used to 
    weight the observations from the original learning data when the trees 
    are refitted. If \code{weight = NULL}, the \code{sampler} will be used.}
  \item{applyfun}{a \code{\link[base]{lapply}}-like function. The default is to use
    \code{\link[base]{lapply}} unless \code{cores} is specified in which case
    \code{\link[parallel]{mclapply}} is used (for multicore computations on platforms
    that support these).}
  \item{cores}{integer. The number of cores to use in multicore computations
    using \code{\link[parallel]{mclapply}} (see above).}
  \item{savetrees}{logical. If \code{TRUE}, trees based on resampled data sets are returned.}
  \item{\dots}{further arguments passed to \code{sampler}.}
}

\details{
  The function \code{stabletree} assesses the stability of tree learners (i.e.,
  recursive partitioning methods) by refitting the tree to resampled versions
  of the learning data. By default, if \code{data = NULL}, the fitting call is
  extracted by \code{\link[stats]{getCall}} to infer the learning data.
  Subsequently, the \code{sampler} generates \code{B} resampled versions
  of the learning data, the tree is regrown with \code{\link[stats]{update}},
  and (if necessary) coerced by \code{\link[partykit]{as.party}}. For each
  of the resampled trees it is queried and stored which variables are selected
  for splitting and what the selected cutpoints are.

  The resulting object of class \code{"stabletree"} comes with a set of 
  standard methods to generic functions including \code{print}, \code{summary}
  for numerical summaries and \code{plot}, \code{barplot}, and \code{image}
  for graphical representations. See \code{\link{plot.stabletree}} for more
  details.
}

\value{
  \code{stabletree} returns an object of class \code{"stabletree"} which is a list with
  the following components:
  \item{call}{the call from the model object \code{x},}
  \item{B}{the number of resampled datasets,}
  \item{sampler}{the \code{sampler} function,}
  \item{vs0}{FIXME better name?,} 
  \item{br0}{FIXME better name?,}
  \item{vs}{FIXME better name?,}
  \item{br}{FIXME better name?,}
  \item{classes}{character vector indicating the classes of all partitioning variables,}
  \item{trees}{a list of tree objects of class \code{"party"}, or \code{NULL}.}
}

\references{
Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.

Philipp M, Zeileis A, Strobl C (2016).
  \dQuote{A Toolkit for Stability Assessment of Tree-Based Learners}. 
  In A. Colubi, A. Blanco, and C. Gatu (Eds.), Proceedings of COMPSTAT 2016 --
  22nd International Conference on Computational Statistics (pp. 315--325). 
  The International Statistical Institute/International Association for
  Statistical Computing. Preprint available at
  \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2016-11}
}

\seealso{\code{\link{plot.stabletree}}, \code{\link[partykit]{as.party}}}

\examples{

\donttest{
## build a simple tree
library("partykit")
m <- ctree(Species ~ ., data = iris)
plot(m)

## investigate stability
set.seed(0)
s <- stabletree(m, B = 500)
print(s)

## variable selection statistics
summary(s)

## show variable selection proportions
barplot(s)

## illustrate variable selections of replications
image(s)

## graphical cutpoint analysis
plot(s)
}

}

\keyword{regression}
