\name{utility.aggregate.mix}

\alias{utility.aggregate.mix}

\title{Mixed aggregation of values and utilities}

\description{
Function to perform a mixed aggregation of values and utilities.
The mixture consists of a weighted mean of the additive, minimum
and Cobb-Douglas aggregation techniques.
}

\usage{
utility.aggregate.mix(u, par)
}

\arguments{
  \item{u}{
numeric vector of values or utilities to be aggregated.
}
  \item{par}{
numeric vector of weights for calculating the weighted mean of the
values provided in the argument \code{u} followed by the three 
weights of the additive, minimum and Cobb-Douglas aggregation 
techniques.
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
Constructing, evaluating and visualizing value and utility functions for decision support, submitted. \cr\cr
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
\code{\link{utility.aggregate.max}}, \cr
\code{\link{utility.aggregate.mult}}.
}

\examples{
utility.aggregate.mix(c(0.2,0.8),par=c(1,1 , 1,0,0))
utility.aggregate.mix(c(0.2,0.8),par=c(1,1 , 0,1,0))
utility.aggregate.mix(c(0.2,0.8),par=c(1,1 , 0,0,1))
utility.aggregate.mix(c(0.2,0.8),par=c(1,1 , 1,1,1))
}

