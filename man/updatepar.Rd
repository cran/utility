\name{updatepar}

\alias{updatepar}

\title{Update Parameters in Node Definitions}

\description{
Generic function to update parameters in all node defintions of the hierarchy defined by the given node.
}

\usage{
updatepar(x, ...)
}

\arguments{
  \item{x}{
node to be updated.
}
  \item{...}{
parameter values can be provided by an additional argument \code{par}.
}
}

\value{
The node or node hierarchy with updated parameters is returned.
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
See \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.cond.create}}, \cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, \cr
\code{\link{utility.conversion.parfun.create}} \cr\cr
for how to construct the nodes and \cr\cr
\code{\link{updatepar.utility.endnode.discrete}} \cr
\code{\link{updatepar.utility.endnode.intpol1d}} \cr
\code{\link{updatepar.utility.endnode.parfun1d}} \cr
\code{\link{updatepar.utility.endnode.intpol2d}} \cr
\code{\link{updatepar.utility.endnode.cond}} \cr
\code{\link{updatepar.utility.aggregation}} \cr
\code{\link{updatepar.utility.conversion.intpol}} \cr
\code{\link{updatepar.utility.conversion.parfun}} \cr
for the updates of the specific nodes.
}


