\name{utility.endnode.intpol2d.create}

\alias{utility.endnode.intpol2d.create}

\title{Construct a two-attribute interpolation end node}

\description{Function to construct a two-attribute interpolation end node.}

\usage{
utility.endnode.intpol2d.create(name.node, 
                                name.attrib, 
                                ranges, 
                                isolines, 
                                u, 
                                names.u      = rep(NA, length(u)), 
                                lead         = 0, 
                                utility      = TRUE, 
                                required     = FALSE, 
                                col          = "black", 
                                shift.levels = 0)
}

\arguments{
  \item{name.node}{
name of the node to be constructed as a character string.
}
  \item{name.attrib}{
names of the attributes on which the value or utility function depends as a vector of two character strings.
}
  \item{ranges}{
list of two numeric vectors with two components each specifying the minimum and the maximum of the range of the corresponding attribute.
}
  \item{isolines}{
list of isoline definitions.
Each definition consists of a list with elements \code{x} and \code{y} that each represents a numeric vector of x- (=first attribute) and y- (second attribute) values to characterize the shape of the isoline.
}
  \item{u}{
numeric vector of the same length as the outer list of the argument \code{isolines} specifying the corresponding values or utilities.
}
  \item{names.u}{
(optional) vector of character strings with names of the components of the numeric vector \code{u} specified above.
Only required to provide acces through a named parameter vector.
}
  \item{lead}{
numeric value specifying which variable is the lead variable for interpolation.
1 indicates linear interpolation between isolines along lines with constant value of the first attribute, 2 along lines with constant values of the second attribute, and zero indicates to take the average of these two interpolation schemes.
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
The function returns the created object of type \code{utility.endnode.intpol2d} with the properties specified in the arguments of the function.
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
\code{\link{print.utility.endnode.intpol2d}}, \cr
\code{\link{summary.utility.endnode.intpol2d}}, \cr
\code{\link{evaluate.utility.endnode.intpol2d}} and \cr
\code{\link{plot.utility.endnode.intpol2d}}. \cr\cr
Create other end nodes with \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.cond.create}}, or \cr
\code{\link{utility.endnode.firstavail.create}}. \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

\examples{
  riparzone_width <-
    utility.endnode.intpol2d.create(
      name.node   = "riparian zone width",
      name.attrib = c("riverbedwidth_m","riparianzonewidth_m"),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,2,10,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead         = 1,
      utility      = FALSE)

  print(riparzone_width)
  plot(riparzone_width)
}
