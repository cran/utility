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

\details{
The aggregation function is defined by
\deqn{u = 1 - \prod_{i=1}^n (1-u_i+\delta)^{w_i} + \delta}{u = 1 - prod (1-ui+delta)^wi + delta}
where \eqn{\delta}{delta} is the last parameter appended to the weights.\cr
The following figure shows examples of the behaviour of this aggregation function and its special case \code{\link{utility.aggregate.revgeo}} for the two-dimensional case:\cr
\if{html}{\figure{aggregationrevgeo.png}{options: width=80\%}}
\if{latex}{\figure{aggregationrevgeo.pdf}{options: width=5in}}
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
utility.aggregate.revgeooff(c(0.2,0.8), par=c(1,1,0.1))
}

