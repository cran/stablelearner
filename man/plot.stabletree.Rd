\name{plot.stabletree}
\alias{plot.stabletree}
\alias{barplot.stabletree}
\alias{image.stabletree}

\title{Visualizing Tree Stability Assessments}

\description{
  Visualizations of tree stability assessments carried out
  via \code{\link{stabletree}}.
}

\usage{
\method{plot}{stabletree}(x, select = order(colMeans(x$vs), decreasing = TRUE), 
  type.breaks = "levels", col.breaks = "red", lty.breaks = "dashed",
  cex.breaks = 0.7,  col.main = c("black", "gray50"), main.uline = TRUE,
  args.numeric = NULL, args.factor = NULL, args.ordered = NULL, main = NULL,
  original = TRUE, \dots)

\method{barplot}{stabletree}(height, main = "Variable selection frequencies",
  xlab = "", ylab = "", horiz = FALSE, col = gray.colors(2),
  names.arg = NULL, names.uline = TRUE, names.diag = TRUE,
  cex.names = 0.9, 
  ylim = if (horiz) NULL else c(0, 100), xlim = if (horiz) c(0, 100) else NULL,
  original = TRUE, \dots)

\method{image}{stabletree}(x, main = "Variable selections",
  ylab = "Repetitions", xlab = "", col = gray.colors(2),
  names.arg = NULL, names.uline = TRUE, names.diag = TRUE, 
  cex.names = 0.9, xaxs = "i", yaxs = "i",
  col.tree = 2, lty.tree = 2, xlim = c(0, length(x$vs0)), ylim = c(0, x$B),
  original = TRUE, \dots)
}

\arguments{
  \item{x, height}{an object of class \code{\link{stabletree}}.}
  \item{original}{logical. Should the original tree information be 
  highlighted?}
  \item{select}{An vector of integer or character values representing the 
  number(s) or the name(s) of the variable(s) to be plotted. By default all 
  variables are plotted. The numbers correspond to the ordering of all 
  partitioning variables used in the call of the fitted model object that was 
  passed to \code{\link{stabletree}}.}
  \item{type.breaks}{A character specifying the type of information added to 
  the lines that represent the splits in the complete data tree. 
  \code{type.breaks = "levels"} adds the level of the splits.
  \code{type.breaks = "nodeids"} adds the nodeid of the splits.
  \code{type.breaks = "breaks"} adds the cutpoint of the splits.
  \code{type.breaks = "none"} suppresses any labeling.}
  \item{col.breaks}{Coloring of the lines and the texts that represent the 
  splits in the complete data tree.}
  \item{lty.breaks}{Type of the lines that represent the splits in the complete 
  data tree.}
  \item{cex.breaks}{Size of the texts that represent the splits in the complete 
  data tree.}
  \item{col.main}{A vector of colors of length two. The first color is used for 
  titles of variables that are selected in the complete data tree. The second 
  color is used for titles of variables that are not selected in the complete 
  data tree.}
  \item{main.uline}{A logical value. If \code{TRUE}, variables selected in the 
  complete data tree are underlined. If \code{FALSE}, underlining is 
  suppressed.}
  \item{args.numeric}{A list of arguments passed to the internal function 
  that is used for plotting a histogram of the cutpoints in numerical splits.
  \code{breaks} is passed to \code{hist} (see \code{\link{hist}} for details). 
  Further arguments in the list are passed to \code{\link{plot.histogram}} such 
  as \code{col}, \code{border} etc.}
  \item{args.factor}{A list of arguments passed to the internal function that 
  is used for plotting an image plot of the cutpoints in categorical splits.
  \code{col} (a vector of two colors) is used to illustrate 
  categorical splits. Please note that the default colors are optimized for 
  color vision deficiency. \code{col.na} defines the color used to highlight 
  missing categories. Further arguments in the list are passed to the function 
  \code{\link{plot.default}}.}
  \item{args.ordered}{A list of arguments passed to the internal function that 
  is used for plotting a barplot of the cutpoints in ordered categorical splits. 
  All arguments in the list are passed to the function 
  \code{\link{barplot.default}}, such as \code{col}, \code{border} etc.}
  \item{\dots}{further arguments passed to plotting functions,
    especially for labeling and annotation.}
  \item{main, xlab, ylab}{character. Annotations of axes and main
    title, respectively.}
  \item{horiz}{A logical value. If \code{FALSE}, the bars are drawn vertically 
  with the first bar to the left. If \code{TRUE}, bars are drawn horizontally 
  with the first at the bottom.}
  \item{col}{A vector of colors of length two used for coloring in the
  \code{barplot} and \code{image} methods. The first color represents selected 
  and the second color represents not selected partitioning variables.}
  \item{names.arg}{A vector of labels to be plotted below each bar (in case of 
  \code{barplot}) or each column (in case of \code{image}). If the argument is omitted, then 
  the variable names are taken as labels.}
  \item{names.uline}{A logical value. If \code{TRUE}, the labels representing 
  variables that are used for splitting in the complete data tree, are 
  underlined. If \code{FALSE}, labels are not underlined.}
  \item{names.diag}{A logical value (omitted if \code{horiz = TRUE}). If 
  \code{TRUE}, the labels are drawn diagonally. If \code{FALSE}, labels are 
  drawn horizontally.}
  \item{cex.names}{Expansion factor for labels.}
  \item{xlim, ylim}{The limits of the plot.}
  \item{xaxs, yaxs}{The style of axis interval calculation to be used (see 
    \code{\link{par}} for details).}
  \item{col.tree, lty.tree}{color and line type to indicate differences from
    the original tree that was resampled.}
}

\details{
  \itemize{
    \item{\code{plot} visualizes the variability of the cutpoints.}
    \item{\code{barplot} visualizes the variable selection frequency.}
    \item{\code{image} visualizes the combinations of variables selected.}
  }
}

\seealso{\code{\link{stabletree}}}

\examples{

\donttest{

## build a tree
library("partykit")
m <- ctree(Species ~ ., data = iris)
plot(m)

## investigate stability
set.seed(0)
s <- stabletree(m, B = 500)

## show variable selection proportions
## with different labels and different ordering
barplot(s)
barplot(s, cex.names = 0.8)
barplot(s, names.diag = FALSE)
barplot(s, names.arg = c("a", "b", "c", "d"))
barplot(s, names.uline = FALSE)
barplot(s, col = c("lightgreen", "darkred"))
barplot(s, horiz = TRUE)

## illustrate variable selections of replications
## with different labels and different ordering
image(s)
image(s, cex.names = 0.8)
image(s, names.diag = FALSE)
image(s, names.arg = c("a", "b", "c", "d"))
image(s, names.uline = FALSE)
image(s, col = c("lightgreen", "darkred"))
image(s, order.row = FALSE)
image(s, order.col = FALSE)

## graphical cutpoint analysis, selecting variable by number and name
## with different numerical of break points
plot(s)
plot(s, select = 3)
plot(s, select = "Petal.Width")
plot(s, args.numeric = list(breaks = 40))

# change labels of splits in complete data tree
plot(s, select = 3, type.labels = "levels")
plot(s, select = 3, type.labels = "nodeids")
plot(s, select = 3, type.labels = "breaks")
plot(s, select = 3, type.labels = "none")

}

}

\keyword{regression}
