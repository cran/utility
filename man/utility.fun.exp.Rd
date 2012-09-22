\name{utility.fun.exp}

\alias{utility.fun.exp}

\title{Exponential function for value or utility functions}

\description{Exponential function for value or utility functions.}

\usage{
utility.fun.exp(attrib, par)
}

\arguments{
  \item{attrib}{
vector of attribute levels to calculate corresponding value or utility.
}
  \item{par}{
Vector of parameters: \cr
\code{par[1]}: absolute risk aversion \cr
\code{par[2]}: minimum of attribute range (default = 0) \cr
\code{par[3]}: maximum of attribute range (default = 1) 
}
}

\details{
The function evaluates the expression \cr
\code{(1-exp(-par[1]*(a-par[2])/(par[3]-par[2])))/(1-exp(-par[1]))}.
}

\value{
Vector of values or utilities corresponding to the attributes passed by argument \code{a}
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
See the node constructors \cr
\code{\link{utility.endnode.intpol1d.create}} and \code{\link{utility.conversion.intpol.create}} \cr
in which this function can be used.
}

\examples{
utility.fun.exp(0:10/10,par=c(2,0,1))
}

