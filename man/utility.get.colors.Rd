\name{utility.get.colors}

\alias{utility.get.colors}

\title{Get Color Corresponding to Specified Value Levels}

\description{Function to get the colors from a given color scheme at specific value levels.}

\usage{
utility.get.colors(u,col=utility.calc.colors())
}

\arguments{
  \item{u}{
value level representing the evaluation ot a value function (this value level has to be between zero and unity).
}
  \item{col}{
color scheme (vector of colors to be used for a division of the interval between zero and unity into equal intervals)).
}
}

\value{
Character vector of colors.
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
See \cr
\code{\link{utility.calc.colors}}
}

\examples{
utility.get.colors(c(0,0.5,1))
}

