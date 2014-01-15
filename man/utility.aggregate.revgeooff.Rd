\name{utility.aggregate.revgeooff}

\alias{utility.aggregate.revgeooff}

\title{Reverse geometric aggregation of values or utilities with offset}

\description{
Function to perform a reverse geometric aggregation 
(unity minus the weighted geometric mean of unity minus the arguments) 
of values or utilities with offset.
}

\usage{
utility.aggregate.revgeooff(u, par)
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
numeric vector of weights for calculating the reverse weighted geometric mean 
of the values provided in the argument \code{u}.
The weights need not be normalized, they will be normalized before use.
In case of missing values in the vector \code{u}, the weights of the
non-missing components will be rescaled to sum to unity.
}
}

\value{
The function returns the aggregated value or utility.
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
Constructor of aggregation node: \cr\cr
\code{\link{utility.aggregation.create}} \cr\cr
Alternative aggregation techniques: \cr\cr
\code{\link{utility.aggregate.add}}, \cr
\code{\link{utility.aggregate.min}}, \cr
\code{\link{utility.aggregate.max}}, \cr
\code{\link{utility.aggregate.cobbdouglas}}, \cr
\code{\link{utility.aggregate.geo}}, \cr
\code{\link{utility.aggregate.geooff}}, \cr
\code{\link{utility.aggregate.revgeo}}, \cr
\code{\link{utility.aggregate.harmo}}, \cr
\code{\link{utility.aggregate.harmooff}}, \cr
\code{\link{utility.aggregate.revharmo}}, \cr
\code{\link{utility.aggregate.revharmooff}}, \cr
\code{\link{utility.aggregate.mult}}, \cr
\code{\link{utility.aggregate.mix}}.
}

\examples{
utility.aggregate.revgeooff(c(0.2,0.8), par=c(1,1))
}

