\name{utility.aggregate.bonusmalus}

\alias{utility.aggregate.bonusmalus}

\title{Bonus-malus aggregation of values or utilities}

\description{
Function to perform an aggregation of valus or utilities that considers some of the 
inputs only as bonus (only considered if value is larger then the aggregated value
of the non bonus or malus input) or malus (only considered if value is smaller then 
the aggregated value of the non bonus or malus input).
}

\usage{
utility.aggregate.bonusmalus(u,par,def.agg="utility.aggregate.add")
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
numeric vector combining the parameters of the default aggregation technique
(see argument \code{def.agg}) with those specifying the bonus-malus behaviour.
The arguments of \code{def.agg}) must match the number of arguments of this
function for the number of inputs reduced to those that are not treated as 
bonus or malus.
This parameter vector is then appended by the parameters characterizing the
bonus-malus behavior. 
This is a parameter vector of the same length as the number of sub-objectives.
Its elements must be NA for the sub-objectives onsidered for the default aggregation
technique, the weights relative to the aggregated value of the non-bonus and non-malus sub-objectives for the sub-objectives to be considered as bonus objectives, and the weights with a negative sign for those to be considered as malus objectives.
Note that the weights of the bonus or malus attributes are relative to the 
aggregated result of the non-bonus and non-malus inputs and the negative signs will only be used for identifying malus sub-objectives and will be eliminated when calculating the weighted mean.
}
  \item{def.agg}{
(optional) character string specifying the name of the function used for aggregation
of the non-bonus and non-malus sub-objectives. 
Note that for use of this aggregation technique in the function 
\code{\link{utility.aggregation.create}}, this argument has to be specified as the
input argument \code{def.agg} (default aggregation) unless it should be additive (default).
}
}

\value{
The function returns the aggregated value or utility.
}

\details{
The aggregation function is defined by
\deqn{u = \frac{\displaystyle u^{\mathrm{agg}}_{i \notin b,i \notin m} + \sum_{\begin{array}{l}i \in b \wedge u_i > u^{\mathrm{agg}}_{i \notin b,i \notin m}\\i \in m \wedge u_i < u^{\mathrm{agg}}_{i \notin b,i \notin m}\end{array}}\mid w_i \mid u_i}{\displaystyle 1 + \sum_{\begin{array}{l}i \in b \wedge u_i > u^{\mathrm{agg}}_{i \notin b,i \notin m}\\i \in m \wedge u_i < u^{\mathrm{agg}}_{i \notin b,i \notin m}\end{array}}\mid w_i \mid}}{u = additive aggregation between the aggregation result of the values or utilities of the non-bonus and non-malus sub-objectives and the values of bonus sub-objectives that are larger than this aggregated value, and the values of malus sub-objectives that are smaller than the aggregated value}
The following figure shows examples of the behaviour of this aggregation function for the two-dimensional case:\cr
\if{html}{\figure{aggregationbonusmalus.png}{options: width=80\%}}
\if{latex}{\figure{aggregationbonusmalus.pdf}{options: width=5in}}
}

\note{
This is the same function as \code{\link{utility.aggregate.cobbdouglas}}
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
utility.aggregate.bonusmalus(c(0.2,0.8), par=c(1,NA,1))
utility.aggregate.bonusmalus(c(0.2,0.8), par=c(1,1,NA))
utility.aggregate.bonusmalus(c(0.2,0.8), par=c(1,NA,-1))
utility.aggregate.bonusmalus(c(0.2,0.8), par=c(1,-1,NA))
}

