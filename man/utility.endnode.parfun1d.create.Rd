\name{utility.endnode.parfun1d.create}

\alias{utility.endnode.parfun1d.create}

\title{Construct a single-attribute parametric function end node}

\description{Function to construct a single-attribute parametric function end node.}

\usage{
utility.endnode.parfun1d.create(name.node, 
                                name.attrib, 
                                range, 
                                name.fun, 
                                par, 
                                names.par    = rep(NA, length(par)), 
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
  \item{name.fun}{
name of the parametric function to be evaluated as a character string.
The parametric function must have the arguments \code{attrib} and \code{par} 
which pass a vector of attribute levels and a vector of parameters
to the function, respectively.
The function has to return a vector of corresponding values or utilities.
}
  \item{par}{
numeric vector of parameter values to be passed to the function specified unter \code{name.fun}.
}
  \item{names.par}{
(optional) vector of parameter names corresponging to the vector of values specified under \code{par}.
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
The function returns the created object of type \code{utility.endnode.parfun1d} with the properties specified in the arguments of the function.
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
\code{\link{print.utility.endnode.parfun1d}}, \cr
\code{\link{summary.utility.endnode.parfun1d}}, \cr
\code{\link{evaluate.utility.endnode.parfun1d}} and \cr
\code{\link{plot.utility.endnode.parfun1d}}. \cr\cr
Create other end nodes with \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.cond.create}}, or \cr
\code{\link{utility.endnode.firstavail.create}}. \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

\examples{
  bedmod_other <-
    utility.endnode.parfun1d.create(
      name.node   = "bed modification other",
      name.attrib = "bedmodfract_percent",
      range       = c(0,100),
      name.fun    = "utility.fun.exp",
      par         = c(-1,100,0),
      required    = FALSE,
      utility     = FALSE)
      
  print(bedmod_other)
  plot(bedmod_other)
}

