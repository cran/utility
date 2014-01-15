\name{evaluate.utility.conversion.intpol}

\alias{evaluate.utility.conversion.intpol}

\title{Evaluate Node and Associated Hierarchy}

\description{
Calculate values or utilities at all nodes of a hierarchy 
for given levels of the attributes.}

\usage{
\method{evaluate}{utility.conversion.intpol}(x, attrib, par = NA, ...)
}

\arguments{
  \item{x}{
node to be evaluated.
}
  \item{attrib}{
numeric vector with labelled components providing the levels 
of a single set of attributes or data frame for which each row
provides such a set of attributes.
}
  \item{par}{
(optional) labelled numeric parameter vector providing parameters to
modify the value or utility function before evaluation.
}
\item{...}{
currently no other arguments are implemented or passed further.
}
}

\value{
Data frame with results of values or utilities 
at all nodes of the hierarchy 
for all provided sets of attribute levels.
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
Peter Reichert <peter.reichert@eawag.ch>
}

\seealso{
\code{\link{utility.conversion.intpol.create}} to create the node, \cr
\code{\link{print.utility.conversion.intpol}} or 
\code{\link{summary.utility.conversion.intpol}} to print its definition, and \cr
\code{\link{plot.utility.conversion.intpol}} to plot the node \cr\cr
and \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.cond.create}}. \cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.parfun.create}} \cr\cr
to create other nodes.
}

\examples{
# see
help(utility)
# for examples.
}

