\name{utility.conversion.parfun.create}

\alias{utility.conversion.parfun.create}

\title{Construct a parametric function conversion node}

\description{Function to construct a node converting values into utilities by a parametric function.}

\usage{
utility.conversion.parfun.create(name.node, 
                                 node, 
                                 name.fun, 
                                 par, 
                                 names.par    = rep(NA, length(par)), 
                                 required     = FALSE, 
                                 col          = "black", 
                                 shift.levels = 0)
}

\arguments{
  \item{name.node}{
name of the node to be constructed as a character string.
}
  \item{node}{
value node that is to be converted into a utility node.
}
  \item{name.fun}{
name of the parametric function to be evaluated as a character string.
The parametric function must have the arguments \code{u} and \code{par} 
which pass a vector of values and a vector of parameters
to the function, respectively.
The function has to return a vector of corresponding utilities.
}
  \item{par}{
numeric vector of parameter values to be passed to the function specified unter \code{name.fun}.
}
  \item{names.par}{
(optional) vector of parameter names corresponging to the vector of values specified under \code{par}.
Only required to provide access to the values through a named parameter vector.
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
The function returns the created object of type \code{utility.conversion.parfun} with the properties specified in the arguments of the function.
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
Print, evaluate and plot the node with \cr\cr
\code{\link{print.utility.conversion.parfun}}, \cr
\code{\link{summary.utility.conversion.parfun}}, \cr
\code{\link{evaluate.utility.conversion.parfun}} and \cr
\code{\link{plot.utility.conversion.parfun}}. \cr\cr
Create other conversion nodes with
\code{\link{utility.conversion.intpol.create}}.
Create end nodes with \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, or \cr
\code{\link{utility.endnode.cond.create}}. \cr\cr
Create aggregation nodes with \cr\cr
\code{\link{utility.aggregation.create}}.
}


