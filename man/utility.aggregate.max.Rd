\name{utility.aggregate.max}

\alias{utility.aggregate.max}

\title{Maximum aggregation of values or utilities}

\description{Function to perform a maximum aggregation of values or utilities.}

\usage{
utility.aggregate.max(u, par = NA)
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
unused argument used for compatibility with other aggregation techniques
that require parameters.
}
}

\value{
maximum of the components of \code{u}.
}

\references{
Short description of the package: \cr\cr
Reichert, P., Schuwirth, N. and Langhans, S., 
Constructing, evaluating and visualizing value and utility functions for decision support, Environmental Modelling & Software, in press. \cr\cr
Textbooks on the use of utility and value functions in decision analysis: \cr\cr
Keeney, R. L. and Raiffa, H. Decisions with Multiple Objectives. Cambridge University Press, 1993. First published in 1976. \cr\cr
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
\code{\link{utility.aggregate.cobbdouglas}}, \cr
\code{\link{utility.aggregate.mult}}, \cr
\code{\link{utility.aggregate.mix}}.
}

\examples{
utility.aggregate.max(c(0.2,0.8))
}

