\name{plot.utility.conversion.intpol}

\alias{plot.utility.conversion.intpol}

\title{Plot Node Definition or Underlying Objectives Hierarchy}

\description{
Plot node definition or underlying objective hierarchy.
}

\usage{
\method{plot}{utility.conversion.intpol}(x, 
     u           = NA, 
     uref        = NA, 
     par         = NA, 
     type        = c("hierarchy", "table", "node", "nodes"), 
     nodes       = NA, 
     col         = utility.calc.colors(), 
     gridlines   = c(0.2, 0.4, 0.6, 0.8), 
     main        = "", 
     cex.main    = 1, 
     cex.nodes   = 1, 
     cex.attrib  = 1, 
     f.reaches   = 0.2,
     f.nodes     = 0.2,
     with.attrib = TRUE,
     levels      = NA,
     plot.val    = TRUE,
     print.val   = TRUE,
     two.lines   = FALSE,
     ...)
}

\arguments{
  \item{x}{
node to be plotted.
}
  \item{u}{
(optional) vector or data frame with elements or columns labelled
according to the nodes of the hierarchy containing values or utilities.
Typically, this will be the complete output or an output row of
the function \cr \code{\link{evaluate.utility.aggregation}}. \cr
This input is only considered if the argument \code{type} is specified
to be either \code{"hierarchy"} or \code{"table"}.
It is then used to color-code the boxes of the hierarchy representing value nodes or the table.
If \code{u} is a data frame with more than one row and the argument \code{type} 
is equal to \code{"hierarchy"}, then
the median and quantile boxes are plotted for value nodes or the expected utility for utility nodes unless the argument \code{main}
contains as many elements as the number of rows of \code{u}.
In the latter case, separate hierarchies with color-coded boxes for value nodes are produced for all rows of \code{u}.
For \code{type} equals \code{"table"}, this argument can be a list of data frames
to make it possible to plot uncertainty ranges from the samples provided
in the list.
}
  \item{uref}{
(optional) vector or data frame with elements or columns labelled
according to the nodes of the hierarchy containing values or utilities.
Typically, this will be the complete output or an output row of
the function \cr \code{\link{evaluate.utility.aggregation}}. \cr 
This input is only considered if the argument \code{type} is specified
to be \code{"hierarchy"}.
It is then used to color-code the upper part of the boxes of the hierarchy
to allow for a comparison with the results provided by the argument \code{u} which are shown in the lower part of the boxes.
}
  \item{par}{
(optional) labelled numeric parameter vector providing parameters to
modify the value or utility function before plotting the node.
Note that this affects only the node definitions plotted if the argument \code{type}
is specified to be \code{"node"} or \code{"nodes"}.
To color-code hierarchies or tables for different parameter values,
the parameters have to be passed to \cr \code{\link{evaluate.utility.aggregation}} \cr before passing the results of this function to this plotting routine.
}
  \item{type}{
(optional) specifies the type of plot to be produced. \cr Options: 
\code{"hierarchy"}, \code{"table"}, \code{"node"} or \code{"nodes"}. \cr
\code{"hierarchy"}: produces a plot of the objectives hierarchy including
color-coded results for values or utilities if these values are provided
by the arguments \code{u} and/or \code{uref}. \cr
\code{"table"}: produces a table with color-coded results for values or
utilities if these values are provided by the argument \code{u}. \cr
\code{"node"}: produces a plot of the definition of the current node. \cr
\code{"nodes"}: produces plots of node definitions for all nodes defined
by the attribute \code{nodes}.
}
  \item{nodes}{
(optional) character vector specifying the nodes for which the definitions
will be plotted or which will be considered in a table. 
The default value of NA indicates that all nodes will be plotted.
This argument only affects the output if the argument \code{type} was indicated to be either \code{"table"} or \code{"nodes"}.
}
  \item{col}{
(optional) character vector of colors to be used to color the interval
between zero and unity in equidistant sections (use repetitions of the
same color if you want to have a non-equidistant color-coding).
This attribute is only used for value nodes and if values are provided
by the arguments \code{u} and/or \code{uref}.
}
  \item{gridlines}{
(optional) numeric vector of levels at which gridlines are plotted
in node definitions.
This attribute is only used if the argument \code{type} is specified to be either \code{"node"} or \code{"nodes"}.
}
  \item{main}{
(optional) title(s) of the plot.
If the argument \code{type} is equal to \code{"hierarchy"} and the a vector of titles with the same length as the number of rows of the argument \code{u} is provided, a color-coded hierarchy is plotted for each row of \code{u}.
Otherwise, the medians and colored boxes indicating 90\% credibility or occurrence ranges are plotted at all nodes.
}
  \item{cex.main}{
(optional) scaling factor for title of the plot.
}
  \item{cex.nodes}{
(optional) scaling factor for node labels used in the plot.
}
  \item{cex.attrib}{
(optional) scaling factor for attribute labels used in the plot.
}
  \item{f.reaches}{
(optional) fraction of the width of the plot reserved for the row labels of the table if the argument \code{type} is equal to \code{"table"}.
}
  \item{f.nodes}{
(optional) fraction of the height of the plot reserved for the column labels of the table if the argument \code{type} is equal to \code{"table"}.
}
  \item{with.attrib}{
(optional) indicates if attributes should be listed if the argument \code{type} is equal to \code{"hierarchy"}.
}
  \item{levels}{
(optional) how many levels of the hierarchy should be plotted (NA means to plot all levels).
}
  \item{plot.val}{
(optional) plot value as a vertical line within the box.
}
  \item{print.val}{
(optional) print value as a number when plotting a table of boxes.
}
  \item{two.lines}{
(optional) choose whether two lines should be used for the labels in the hierarchy plot.
}
  \item{\dots}{
additional arguments passed to the R plotting routine.
}
}

\note{
Note that the plotting routines \cr
\code{\link{plot.utility.conversion.parfun}} \cr
\code{\link{plot.utility.aggregation}} \cr
are exactly the same so that all hierarchies can be plotted with exactly
the same commands irrespective of the type of the top-level node.
}

\references{
Short description of the package: \cr\cr
Reichert, P., Schuwirth, N. and Langhans, S., 
Constructing, evaluating and visualizing value and utility functions for decision support, Environmental Modelling & Software 46, 283-291, 2013. \cr\cr
Textbooks on the use of utility and value functions in decision analysis: \cr\cr
Keeney, R. L. and Raiffa, H. Decisions with Multiple Objectives - Preferences and Value Tradeoffs. John Wiley & Sons, 1976. \cr\cr
Eisenfuehr, F., Weber, M. and Langer, T., Rational Decision Making, Springer, Berlin, 2010.
}

\author{
Peter Reichert <peter.reichert@emeriti.eawag.ch>
}

\seealso{
See \code{\link{utility.conversion.intpol.create}} for how to construct such a node and \code{\link{evaluate.utility.conversion.intpol}} for how to evaluate the node. \cr\cr
See \code{\link{utility.calc.colors}} for an example of how to construct color schemes and \code{\link{utility.get.colors}} for how to get colors for specifed value levels.
}

\examples{
# see
help(utility)
# for examples.
}

