\name{evaluate}

\alias{evaluate}

\title{Evaluate Node and Associated Hierarchy}

\description{
Generic function to calculate values or utilities at all nodes of a hierarchy 
for given levels of the attributes.
}

\usage{
evaluate(x, ...)
}

\arguments{
  \item{x}{
node to be evaluated.
}
\item{...}{
attribute levels have to be provided as an additional argument \code{attrib};
parameter values can optionally be provided as an additional argument \code{par}.
}
}

\value{
Data frame with results of values or utilities 
at all nodes of the hierarchy 
for all provided sets of attribute levels.
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
See \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.cond.create}}, \cr
\code{\link{utility.aggregation.create}}, \cr
\code{\link{utility.conversion.intpol.create}}, \cr
\code{\link{utility.conversion.parfun.create}} \cr\cr
to create the nodes to be evaluated.
}

\examples{
# see
help(utility)
# for examples.
}

