\name{utility.conversion.intpol.create}

\alias{utility.conversion.intpol.create}

\title{Construct an interpolation conversion node}

\description{Function to construct a node converting values into utilities by interpolation.}

\usage{
utility.conversion.intpol.create(name.node, 
                                 node, 
                                 x, 
                                 u, 
                                 names.x      = rep(NA, length(x)), 
                                 names.u      = rep(NA, length(u)), 
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
  \item{x}{
numeric vector of values for which the utility is known.
}
  \item{u}{
numeric vector of utilities corresponding to the values given in the previous argument \code{x}.
}
  \item{names.x}{
(optional) vector of character strings with names of the components of the numeric vector \code{x} specified above.
Only required to provide access to the values through a named parameter vector.
}
  \item{names.u}{
(optional) vector of character strings with names of the components of the numeric vector \code{u} specified above.
Only required to provide acces through a named parameter vector.
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
The function returns the created object of type \code{utility.conversion.intpol1} with the properties specified in the arguments of the function.
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
Print, evaluate and plot the node with \cr\cr
\code{\link{print.utility.conversion.intpol}}, \cr
\code{\link{summary.utility.conversion.intpol}}, \cr
\code{\link{evaluate.utility.conversion.intpol}} and \cr
\code{\link{plot.utility.conversion.intpol}}. \cr\cr
Create other conversion nodes with
\code{\link{utility.conversion.parfun.create}}.
Create end nodes with \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, or \cr
\code{\link{utility.endnode.cond.create}}. \cr\cr
Create aggregation nodes with \cr\cr
\code{\link{utility.aggregation.create}}.
}


