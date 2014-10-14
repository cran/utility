\name{utility.endnode.cond.create}

\alias{utility.endnode.cond.create}

\title{Construct a conditional end node}

\description{
Function to construct a node that makes a choice between given end 
nodes based on the levels of discrete attributes.
}

\usage{
utility.endnode.cond.create(name.node, 
                            attrib.levels, 
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
  \item{attrib.levels}{
data frame with attribute names as column names and all discrete attribute level combinations in the rows.
This may be a dependence on any number of attributes.
As combinatorics can lead to a very large number of possible combinations, the node should not depend on a too large number of attributes, in particular if each attribute has many different levels expressed by numbers or character strings.
}
  \item{nodes}{
list of the length of the number of columns of the data frame specifed as argument \code{attrib.levels} above containing the nodes to be associated with the attribute level combinations specified in the rows of \code{attrib.levels}.
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
The function returns the created object of type \code{utility.endnode.cond} with the properties specified in the arguments of the function.
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
\code{\link{print.utility.endnode.cond}}, \cr
\code{\link{summary.utility.endnode.cond}}, \cr
\code{\link{evaluate.utility.endnode.cond}} and \cr
\code{\link{plot.utility.endnode.cond}}. \cr\cr
Create other end nodes with  \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, or \cr
\code{\link{utility.endnode.firstavail.create}}. \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

\examples{
  bedmod_riprap <-
    utility.endnode.intpol1d.create(
      name.node   = "bed modification riprap",
      name.attrib = "bedmodfract_percent",
      range       = c(0,100),
      x           = c(0,10,30,100),
      u           = c(1,0.775,0.5625,0.24),
      required    = FALSE,
      utility     = FALSE)
  
  bedmod_other <-
    utility.endnode.intpol1d.create(
      name.node   = "bed modification other",
      name.attrib = "bedmodfract_percent",
      range       = c(0,100),
      x           = c(0,10,30,100),
      u           = c(1,0.775,0.5625,0),
      required    = FALSE,
      utility     = FALSE)
  
  bedmod <-
    utility.endnode.cond.create(
      name.node     = "bed modification",
      attrib.levels = data.frame(bedmodtype_class=
                                 c("riprap","other")),
      nodes         = list(bedmod_riprap,bedmod_other),
      required      = FALSE,
      utility       = FALSE)

  print(bedmod)
  plot(bedmod)
}
