\name{utility.calc.colors}

\alias{utility.calc.colors}

\title{Color Scheme for Value Functions}

\description{Function to calculate a color scheme for value functions.}

\usage{
utility.calc.colors(n = 5)
}

\arguments{
  \item{n}{
number of colors.
}
}

\details{
For \code{n = 5} this function produces the standard colors red, orange, yellow, green and blue as used in river assessment programs.
These colors are provided in a lighter version to improve readability of black text in front of the colored background.
For large values of \code{n} quasi-continuous transitions are defined between these colors.
Any other vector of colors can be used by the plotting routines.
}

\value{
Character vector of colors.
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
See \cr
\code{\link{plot.utility.endnode.discrete}} \cr
\code{\link{plot.utility.endnode.intpol1d}} \cr
\code{\link{plot.utility.endnode.parfun1d}} \cr
\code{\link{plot.utility.endnode.intpol2d}} \cr
\code{\link{plot.utility.endnode.cond}} \cr
\code{\link{plot.utility.aggregation}} \cr
\code{\link{plot.utility.conversion.intpol}} \cr
\code{\link{plot.utility.conversion.parfun}} \cr
for the use of such color vectors in plotting functions and \cr
\code{\link{utility.get.colors}} \cr
for getting colors corresponding to specified values.
}

\examples{
utility.calc.colors(5)
utility.calc.colors(100)
}

