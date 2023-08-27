\name{utility.aggregate.mult}

\alias{utility.aggregate.mult}

\title{Multiplicative aggregation of values or utilities}

\description{Function to perform a multiplicative aggregation of values or utilities.}

\usage{
utility.aggregate.mult(u, par)
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
numeric vector of weights for calculating the multiplicative
combination of the values provided in the argument \code{u}.
Note that for this aggregation technique, the result depends on the sum of the weights that need not be unity.
}
}

\value{
numeric value corresponding to the multiplicative aggregation
of the values provided in the vector \code{u}.
}

\details{
To derive the aggregated value, we first solve the implicit equation
\deqn{k + 1 = \prod_{i=1}^n (1 + k w_i)}{k + 1 = prod ( 1 + k wi )}
for \eqn{k} and then calculate the aggregated value as
\deqn{u = \frac{\displaystyle \prod_{i=1}^n (1 + k w_i u_i)}{k}}{u = (prod(1 + k ki  ui) - 1) / k}
See Keeney and Raiffa, Decisions with multiple objectives, 1976, pp. 307, 347-348 for details.\cr
The following figure shows examples of the behaviour of this aggregation function for the two-dimensional case:\cr
\if{html}{\figure{aggregationmult.png}{options: width=80\%}}
\if{latex}{\figure{aggregationmult.pdf}{options: width=5in}}
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
Peter Reichert <peter.reichert@emeriti.eawag.ch>
}

\seealso{
Constructor of aggregation node: \cr\cr
\code{\link{utility.aggregation.create}} \cr\cr
Aggregation techniques provided by uncsim: \cr\cr
\code{\link{utility.aggregate.add}} for additive aggregation (weighted arithmetic mean), \cr
\code{\link{utility.aggregate.min}} for minimum aggregation, \cr
\code{\link{utility.aggregate.max}} for maximum aggregation, \cr
\code{\link{utility.aggregate.geo}} or \code{\link{utility.aggregate.cobbdouglas}} for geometric or Cobb-Douglas aggregation (weighted geometric mean), \cr
\code{\link{utility.aggregate.geooff}} for geometric aggregation with offset, \cr
\code{\link{utility.aggregate.revgeo}} for reverse geometric aggregation, \cr
\code{\link{utility.aggregate.revgeooff}} for reverse geometric aggregation with offset, \cr
\code{\link{utility.aggregate.harmo}} for harmonic aggregation (weighted harmonic mean), \cr
\code{\link{utility.aggregate.harmooff}} for harmonic aggregation with offset, \cr
\code{\link{utility.aggregate.revharmo}} for reverse harmonic aggregation, \cr
\code{\link{utility.aggregate.revharmooff}} for reverse harmonic aggregation with offset, \cr
\code{\link{utility.aggregate.mult}} for multiplicative aggregation, \cr
\code{\link{utility.aggregate.mix}} for a mixture of additive, minimum, and geometric aggregation, \cr
\code{\link{utility.aggregate.addmin}} for a mixture of additive and minimum aggregation. \cr
\code{\link{utility.aggregate.addpower}} for additive power aggregation (weighted power mean), \cr
\code{\link{utility.aggregate.revaddpower}} for reverse additive power aggregation, \cr
\code{\link{utility.aggregate.addsplitpower}} for splitted additive power aggregation, \cr
\code{\link{utility.aggregate.revaddsplitpower}} for reverse splitted additive power aggregation, \cr
\code{\link{utility.aggregate.bonusmalus}} for an aggregation technique that considers some of the values or utilities of sub-objectives only as bonus or malus. \cr
}

\examples{
utility.aggregate.mult(c(0.2,0.8),par=c(0.3,0.3))
}

