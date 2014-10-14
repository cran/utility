\name{utility.aggregate.addmin}

\alias{utility.aggregate.addmin}

\title{Mixture of additive and minimum aggregation}

\description{
Function to perform a mixture of additive and minimum aggregation.
The parameter vector must contain the weights for additive aggregation 
followed by the weight of additive aggregation.
The weight for minimum aggregation is then unity minus the weight for 
additive aggregation.
If this additional weight is zero, we return to minimum aggregation,
if it is unity, we will have additive aggregation.
}

\usage{
utility.aggregate.addmin(u, par)
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
numeric vector of weights for additive aggregation appended by the weight for additive aggregation.
The weight for minimum aggregation is then unity minus the weight for 
additive aggregation.
If this additional weight is zero, we return to minimum aggregation,
if it is unity, we will have additive aggregation.
The weights for additive aggregation need not be normalized, they will be 
normalized before use.
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
Description of aggregation techniques: \cr\cr
Langhans, S.D., Reichert, P. and Schuwirth, N.,
The method matters: A guide for indicator aggregation in ecological assessments. 
Ecological Indicators 45, 494-507, 2014. \cr\cr
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
\code{\link{utility.aggregate.revgeooff}}, \cr
\code{\link{utility.aggregate.harmo}}, \cr
\code{\link{utility.aggregate.harmooff}}, \cr
\code{\link{utility.aggregate.revharmo}}, \cr
\code{\link{utility.aggregate.revharmooff}}, \cr
\code{\link{utility.aggregate.mult}}, \cr
\code{\link{utility.aggregate.mix}}.
}

\examples{
utility.aggregate.addmin(c(0.2,0.8), par=c(1,1,0.5))
}

