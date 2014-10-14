\name{summary.utility.aggregation}

\alias{summary.utility.aggregation}

\title{Print Summary of Definitions of Node and Associated Hierarchy}

\description{
Print summary of definition of node and associated hierarchy.
}

\usage{
\method{summary}{utility.aggregation}(object, ...)
}

\arguments{
  \item{object}{
node of which a summary is to be printed.
}
  \item{...}{
currently no other arguments are implemented or passed further.
}
}

\note{
In the current version of the package, the methods \code{print} and \code{summary}
provide the same output.
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
See \code{\link{utility.aggregation.create}} for how to construct such a node.
}

\examples{
# see
help(utility)
# for examples.
}
