\name{utility.endnode.firstavail.create}

\alias{utility.endnode.firstavail.create}

\title{Construct an end node to get the results of the first available sub-node}

\description{
Function to construct a node that returns the results of the first sub-node for
which results are available.
}

\usage{
utility.endnode.firstavail.create(name.node, 
                                  nodes, 
                                  utility      = TRUE, 
                                  required     = FALSE, 
                                  col          = "black", 
                                  shift.levels = 0)
}

\arguments{
  \item{name.node}{
name of the node to be constructed as a character string.
}
  \item{nodes}{
list of nodes to be tried.
}
  \item{utility}{
(optional) logical variable indicating if a value function (\code{FALSE}) or a utility function (\code{TRUE}) is created.
Default value is \code{TRUE}.
}
  \item{required}{
(optional) logical variable indicating if the value of this node is required for aggregation at the next higher level.
If this variable is \code{TRUE}, aggregation at the next higher level is not possible if this node returns NA.
Default value is \code{FALSE}.
}
  \item{col}{
(optional) color used for plotting the bounding box of the node in the objective hierarchy.
Default value is \code{"black"}.
}
  \item{shift.levels}{
(optional) number of hierarchical levels by which the node in the objective hierarchy is shifted to make a branch fit better to other branches.
Default value is \code{0}.
}
}

\value{
The function returns the created object of type \code{utility.endnode.firstavail} with the properties specified in the arguments of the function.
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
Print, evaluate and plot the node with \cr\cr
\code{\link{print.utility.endnode.firstavail}}, \cr
\code{\link{summary.utility.endnode.firstavail}}, \cr
\code{\link{evaluate.utility.endnode.firstavail}} and \cr
\code{\link{plot.utility.endnode.firstavail}}. \cr\cr
Create other end nodes with  \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, or \cr
\code{\link{utility.endnode.cond.create}}, \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

