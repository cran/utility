\name{utility.endnode.discrete.create}

\alias{utility.endnode.discrete.create}

\title{Construct a discrete value or utility end node}

\description{Function to construct a discrete value or utility end node.}
\usage{
utility.endnode.discrete.create(name.node, 
                                attrib.levels, 
                                u, 
                                names.u      = rep(NA, length(u)), 
                                utility      = TRUE, 
                                required     = FALSE, 
                                col          = "black", 
                                shift.levels = 0)
}

\arguments{
  \item{name.node}{
name of the node to be constructed as a character string.
}
  \item{attrib.levels}{
data frame with attribute names as column names and all discrete attribute level combinations in the rows.
This may be a dependence on any number of attributes.
As combinatorics can lead to a very large number of possible combinations, the node should not depend on a too large number of attributes, in particular if each attribute has many different levels expressed by numbers or character strings.
}
  \item{u}{
numeric vector of the length of the number of columns of the data frame specifed as argument \code{attrib.levels} above specifying the values or utilities corresponding to the rows of \code{attrib.levels}.
}
  \item{names.u}{
(optional) vector of character strings with names of the components of the numeric vector \code{u} specified above.
Only required to provide access to the values through a named parameter vector.
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
The function returns the created object of type \code{utility.endnode.discrete} with the properties specified in the arguments of the function.
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
\code{\link{print.utility.endnode.discrete}}, \cr
\code{\link{summary.utility.endnode.discrete}}, \cr
\code{\link{evaluate.utility.endnode.discrete}} and \cr
\code{\link{plot.utility.endnode.discrete}}. \cr\cr
Create other end nodes with \cr\cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, or \cr
\code{\link{utility.endnode.cond.create}}. \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

\examples{
  widthvar <- 
    utility.endnode.discrete.create(
      name.node     = "width variability",
      attrib.levels = data.frame(widthvariability_class=
                                 c("high","moderate","none")),
      u             = c(1,0.4125,0),
      names.u       = c("u.high","u_moderate","u.none"),
      required      = FALSE,
      utility       = FALSE)

  print(widthvar)
  plot(widthvar)
}
