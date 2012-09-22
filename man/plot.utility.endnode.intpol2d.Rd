\name{plot.utility.endnode.intpol2d}

\alias{plot.utility.endnode.intpol2d}

\title{Plot Node Definition}

\description{
Plot node definition.
}

\usage{
\method{plot}{utility.endnode.intpol2d}(x, 
     par       = NA, 
     col       = utility.calc.colors(utility.num.colors), 
     gridlines = c(0.2, 0.4, 0.6, 0.8), 
     main      = "", 
     cex.main  = 1,
     ...)
}

\arguments{
  \item{x}{
node to be plotted.
}
  \item{par}{
(optional) labelled numeric parameter vector providing parameters to
modify the value or utility function before plotting the node.
}
  \item{col}{
(optional) character vector of colors to be used to color the interval
between zero and unity in equidistant sections (use repetitions of the
same color if you want to have a non-equidistant color-coding).
This attribute is only used for value nodes.
}
  \item{gridlines}{
(optional) numeric vector of levels at which gridlines are plotted
in the node definition.
Not used for this type of node.
}
  \item{main}{
(optional) title of the plot.
}
  \item{cex.main}{
(optional) scaling factor for title of the plot.
}
  \item{\dots}{
additional arguments passed to the R plotting routine.
}
}

\note{
Note that the plotting routines for the other end nodes \cr
\code{\link{plot.utility.endnode.discrete}} \cr
\code{\link{plot.utility.endnode.parfun1d}} \cr
\code{\link{plot.utility.endnode.intpol2d}} \cr
\code{\link{plot.utility.endnode.cond}} \cr
are as far as possible the same so that all end nodes can be plotted with the same commands irrespective of the type of the end node.
}

\references{
Short description of the package: \cr\cr
Reichert, P., Schuwirth, N. and Langhans, S., 
Constructing, evaluating and visualizing value and utility functions for decision support, submitted. \cr\cr
Textbooks on the use of utility and value functions in decision analysis: \cr\cr
Keeney, R. L. and Raiffa, H. Decisions with Multiple Objectives. Cambridge University Press, 1993. First published in 1976. \cr\cr
Eisenfuehr, F., Weber, M. and Langer, T., Rational Decision Making, Springer, Berlin, 2010.
}

\author{
Peter Reichert <peter.reichert@eawag.ch>
}

\seealso{
See \code{\link{utility.endnode.intpol1d.create}} for how to construct such a node and \code{\link{evaluate.utility.endnode.intpol1d}} for how to evaluate the node.
}

\examples{
# see
help(utility)
# for examples.
}
