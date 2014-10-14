\name{updatepar.utility.endnode.firstavail}

\alias{updatepar.utility.endnode.firstavail}

\title{Update Parameters in Node Definitions}

\description{
Update parameters in all node defintions used to define the node.
}

\usage{
\method{updatepar}{utility.endnode.firstavail}(x, par=NA, ...)
}

\arguments{
  \item{x}{
node to be updated.
}
  \item{par}{
parameter vector with labelled parameters to be updated.
}
  \item{...}{
currently no other arguments are implemented or passed further.
}
}

\value{
The node with updated parameters is returned.
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
See \code{\link{utility.endnode.firstavail.create}} for how to construct such a node and \cr\cr
\code{\link{updatepar.utility.endnode.discrete}} \cr
\code{\link{updatepar.utility.endnode.intpol1d}} \cr
\code{\link{updatepar.utility.endnode.parfun1d}} \cr
\code{\link{updatepar.utility.endnode.intpol2d}} \cr
\code{\link{updatepar.utility.endnode.cond}} \cr
\code{\link{updatepar.utility.aggregation}} \cr
\code{\link{updatepar.utility.conversion.intpol}} \cr
\code{\link{updatepar.utility.conversion.parfun}} \cr
for analogous updates of other nodes
}


