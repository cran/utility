\name{utility.aggregate.add}

\alias{utility.aggregate.add}

\title{Additive aggregation of values or utilities}

\description{
Function to perform an additive aggregation (weighted mean) 
of values or utilities.
}

\usage{
utility.aggregate.add(u, par)
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
numeric vector of weights for calculating the weighted mean of the
values provided in the argument \code{u}.
The weights need not be normalized, they will be normalized before use.
In case of missing values in the vector \code{u}, the weights of the
non-missing components will be rescaled to sum to unity.
}
}

\value{
numeric value representing the weighted mean of the components 
of \code{u}.
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
\code{\link{utility.aggregate.min}}, \cr
\code{\link{utility.aggregate.cobbdouglas}}, \cr
\code{\link{utility.aggregate.max}}, \cr
\code{\link{utility.aggregate.mult}}, \cr
\code{\link{utility.aggregate.mix}}.
}

\examples{
utility.aggregate.add(c(0.2,0.8), par=c(1,1))
}

