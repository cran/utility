\name{utility.endnode.classcounts.create}

\alias{utility.endnode.classcounts.create}

\title{Construct an end node that evaluates counts in different quality classes.}

\description{
Function to construct a node that evaluates counts in different quality classes by assigning
a value to the highest class with counts > 0 and optionally increments this value with the
counts in this and lower classses and given increments per count unit.
The user can choose whether these increments can lead to a value higher than the basic value
for counts in the next better class or if this value limits the potential increase.
}

\usage{
utility.endnode.classcounts.create(name.node,        # character(1)
                                   name.attrib,      # character(n)
                                   u.max.inc,        # list (n) of vect (>=1)
                                   names.u.max.inc = list(),
                                   exceed.next     = TRUE,
                                   utility         = TRUE,
                                   required        = FALSE,
                                   col             = "black",
                                   shift.levels    = 0)
}

\arguments{
  \item{name.node}{
name of the node to be constructed as a character string.
}
  \item{name.attrib}{
vector of names of attributes (counts in classes in decreasing order of value of the classes).
}
  \item{u.max.inc}{
list of vectors specifying the basic value and the increments for each class and classes of less value.
}
  \item{names.u.max.inc}{
(optional) list of vectors of names of parameters u.max.inc.
}
  \item{exceed.next}{
(optional) logical variable to indicate whether the level corresponding to the species
classified at the next higher level can be exceeded with increments; 
default value is TRUE.
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
The function returns the created object of type \code{utility.endnode.classcounts} with the properties specified in the arguments of the function.
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
\code{\link{print.utility.endnode.classcounts}}, \cr
\code{\link{summary.utility.endnode.classcounts}}, \cr
\code{\link{evaluate.utility.endnode.classcounts}} and \cr
\code{\link{plot.utility.endnode.classcounts}}. \cr\cr
Create other end nodes with  \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, or \cr
\code{\link{utility.endnode.cond.create}}, \cr
\code{\link{utility.endnode.firstavail.create}}, \cr\cr
Create other types of nodes with \cr\cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

