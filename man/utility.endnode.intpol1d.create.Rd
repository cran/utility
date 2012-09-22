\name{utility.endnode.intpol1d.create}

\alias{utility.endnode.intpol1d.create}

\title{Construct a single-attribute interpolation end node}
\description{Function to construct a single-attribute interpolation end node.}

\usage{
utility.endnode.intpol1d.create(name.node, 
                                name.attrib, 
                                range, 
                                x, 
                                u, 
                                names.x      = rep(NA, length(x)), 
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
  \item{name.attrib}{
name of the attribute on which the value or utility function depends as a character string.
}
  \item{range}{
numeric vector with two components specifying the minimum and the maximum of the attribute range.
}
  \item{x}{
numeric vector of attribute values for which the value or utility is known.
}
  \item{u}{
numeric vector of values or utilities corresponding to the attribute values given in the previous argument \code{x}.
}
  \item{names.x}{
(optional) vector of character strings with names of the components of the numeric vector \code{x} specified above.
Only required to provide access to the values through a named parameter vector.
}
  \item{names.u}{
(optional) vector of character strings with names of the components of the numeric vector \code{u} specified above.
Only required to provide acces through a named parameter vector.
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
The function returns the created object of type \code{utility.endnode.intpol1d} with the properties specified in the arguments of the function.
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
\code{\link{print.utility.endnode.intpol1d}}, \cr
\code{\link{summary.utility.endnode.intpol1d}}, \cr
\code{\link{evaluate.utility.endnode.intpol1d}} and \cr
\code{\link{plot.utility.endnode.intpol1d}}. \cr\cr
Create other end nodes with \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, or \cr
\code{\link{utility.endnode.cond.create}}. \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

\examples{
  bedmod_other <-
    utility.endnode.intpol1d.create(
      name.node   = "bed modification other",
      name.attrib = "bedmodfract_percent",
      range       = c(0,100),
      x           = c(0,10,30,100),
      u           = c(1,0.775,0.5625,0),
      required    = FALSE,
      utility     = FALSE)
      
  print(bedmod_other)
  plot(bedmod_other)
}
