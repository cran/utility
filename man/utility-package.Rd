\name{utility-package}

\alias{utility-package}

\alias{utility}

\docType{package}

\title{Construct, Evaluate and Plot Value and Utility Functions}

\description{Construct and plot objective hierarchies and associated value and utility functions. Evaluate the values and utilities and visualize the results as colored objective hierarchies or tables. Visualize uncertainty by plotting median and quantile intervals within the nodes of the objective hierarchy. Get numerical results of the evaluations in standard \R data types for further processing.}

\details{
\tabular{ll}{
Package: \tab utility\cr
Type: \tab Package\cr
Version: \tab 1.4.5\cr
Date: \tab 2020-03-09\cr
License: \tab GPL-3\cr
}
An objective hierarchy and an associated value or utility function 
is constructed by constructing the nodes of the hierarchy starting 
from the end nodes and proceeding to the higher hierarchies.
Five types of end nodes are distinguished:
End nodes of the class \code{utility.endnode.discrete} define a value 
or utility function for an attribute that has a finite number of
discrete numeric or non-numeric levels.
End nodes of the classes \code{utility.endnode.intpol1d} and 
\code{utility.endnode.parfun1d} implement single-attribute 
value or utility functions that accept a continuous argument.
The first of these functions allows the user to specify 
attribute-value pairs and performs linear interpolation between these
points.
The second function allows the user to specify any parameteric
function that is implemented as a function in R.
End nodes of the class \code{utility.endnode.intpol2d} implement
interpolated value or utility functions that are based on two
attributes.
End nodes of the class \code{utility.endnode.cond} implement
value or utility functions that assign different value or utility
functions to a finite set of attribute combinations.
End nodes of the class \code{utility.endnode.firstavail} implement
value or utility functions that try to evaluate a list of nodes and 
return the value of the first node that could be evaluated based
on the provided attribute data.
Finally, end nodes of the class \code{utility.endnode.classcounts} implement
value or utility functions that value counts e.g. of species of different classes
by assigning a basic value for the occurrence of at least one species of the best class
and incrementing this value by multiplicities of species of this class and of the next lower class.
These end nodes can be implemented by using the following constructors. \cr\cr 
\code{\link{utility.endnode.discrete.create}} \cr
\code{\link{utility.endnode.intpol1d.create}} \cr
\code{\link{utility.endnode.parfun1d.create}} \cr
\code{\link{utility.endnode.intpol2d.create}} \cr
\code{\link{utility.endnode.cond.create}} \cr
\code{\link{utility.endnode.firstavail.create}} \cr
\code{\link{utility.endnode.classcounts.create}} \cr\cr
To advance to higher hierarchical levels, values or utilities at lower
levels must be aggregated to the next higher level.
This is done ab aggregation nodes of the class \code{utility.aggregation}.
Such nodes can be implemented by using the following constructor: \cr\cr
\code{\link{utility.aggregation.create}} \cr\cr
Finally, to provide decision support under uncertainty, values at
an adequate level of the objectives hierarchy must be converted to
utilities by accounting for the risk attitude of the decision maker.
Similar to the single-attribute value or utility functions,
this can either be done by linear interpolation with a node of the class
\code{utility.conversion.intpol} or by using a parametric funciton in
a node of the class \code{utility.conversion.parfun}.
These conversion nodes can be implemented by the constructors: \cr\cr
\code{\link{utility.conversion.intpol.create}} \cr
\code{\link{utility.conversion.parfun.create}} \cr\cr
The definition of the objective hierarchy and the associated value 
and utility function can then be listed or visualized by using the 
generic functions \cr\cr
\code{print} \cr
\code{summary} \cr
\code{plot} \cr\cr 
which automaticall call the implementation corresponding to the node specified
as the first argument: \cr\cr
\code{\link{print.utility.endnode.discrete}} \cr
\code{\link{print.utility.endnode.intpol1d}} \cr
\code{\link{print.utility.endnode.parfun1d}} \cr
\code{\link{print.utility.endnode.intpol2d}} \cr
\code{\link{print.utility.endnode.cond}} \cr
\code{\link{print.utility.endnode.firstavail}} \cr
\code{\link{print.utility.endnode.classcounts}} \cr
\code{\link{print.utility.aggregation}} \cr
\code{\link{print.utility.conversion.intpol}} \cr
\code{\link{print.utility.conversion.parfun}} \cr\cr
\code{\link{summary.utility.endnode.discrete}} \cr
\code{\link{summary.utility.endnode.intpol1d}} \cr
\code{\link{summary.utility.endnode.parfun1d}} \cr
\code{\link{summary.utility.endnode.intpol2d}} \cr
\code{\link{summary.utility.endnode.cond}} \cr
\code{\link{summary.utility.endnode.firstavail}} \cr
\code{\link{summary.utility.endnode.classcounts}} \cr
\code{\link{summary.utility.aggregation}} \cr
\code{\link{summary.utility.conversion.intpol}} \cr
\code{\link{summary.utility.conversion.parfun}} \cr\cr
\code{\link{plot.utility.endnode.discrete}} \cr
\code{\link{plot.utility.endnode.intpol1d}} \cr
\code{\link{plot.utility.endnode.parfun1d}} \cr
\code{\link{plot.utility.endnode.intpol2d}} \cr
\code{\link{plot.utility.endnode.cond}} \cr
\code{\link{plot.utility.endnode.firstavail}} \cr
\code{\link{plot.utility.endnode.classcounts}} \cr
\code{\link{plot.utility.aggregation}} \cr
\code{\link{plot.utility.conversion.intpol}} \cr
\code{\link{plot.utility.conversion.parfun}} \cr\cr
The value or utility function can then be evaluated by applying the 
generic function \cr\cr
\code{evaluate} \cr\cr
that again calls automatically the corresponding class-specific function \cr\cr
\code{\link{evaluate.utility.endnode.discrete}} \cr
\code{\link{evaluate.utility.endnode.intpol1d}} \cr
\code{\link{evaluate.utility.endnode.parfun1d}} \cr
\code{\link{evaluate.utility.endnode.intpol2d}} \cr
\code{\link{evaluate.utility.endnode.cond}} \cr
\code{\link{evaluate.utility.endnode.firstavail}} \cr
\code{\link{evaluate.utility.endnode.classcounts}} \cr
\code{\link{evaluate.utility.aggregation}} \cr
\code{\link{evaluate.utility.conversion.intpol}} \cr
\code{\link{evaluate.utility.conversion.parfun}} \cr\cr
This function requires the provision of observed or predicted attributes
of the valued system and returns the corresponding values or utilities
of all nodes of the hierarchy.
These results can then be visualized by providing them to the generic 
function  \cr\cr
\code{plot} \cr\cr
in addition to the definition of the objective hierarchy stored in the
variable corresponding to the highest node of the hierarchy.
Again, this function automatically calls the correct class-specific
implementation (the root of the hierarchy will be an aggregation 
or a conversion node, not an end node): \cr\cr
\code{\link{plot.utility.aggregation}} \cr
\code{\link{plot.utility.conversion.intpol}} \cr
\code{\link{plot.utility.conversion.parfun}} \cr\cr
This proceedure guarantees easy handling with the simple commands
\code{print}, \code{summary}, \code{evaluate}, and \code{plot} and
the specific function descriptions provided above are only required
to check advanced attributes.
}

\author{
Peter Reichert <peter.reichert@eawag.ch> 
with contributions by Nele Schuwirth <nele.schuwirth@eawag.ch>

Maintainer: Peter Reichert <peter.reichert@eawag.ch>
}

\references{
Short description of the package: \cr\cr
Reichert, P., Schuwirth, N. and Langhans, S., 
Constructing, evaluating and visualizing value and utility functions for decision support, Environmental Modelling & Software 46, 283-291, 2013. \cr\cr
Textbooks on the use of utility and value functions in decision analysis: \cr\cr
Keeney, R. L. and Raiffa, H. Decisions with Multiple Objectives - Preferences and Value Tradeoffs. John Wiley & Sons, 1976. \cr\cr
Eisenfuehr, F., Weber, M. and Langer, T., Rational Decision Making, Springer, Berlin, 2010.
}

\keyword{decision analysis; objectives hierarchy; value function; utility function}

\examples{
# define discrete end node for width variability
# (attribute "widthvariability_class" with levels "high", 
# "moderate" and "none")

widthvar <- 
  utility.endnode.discrete.create(
    name.node     = "width variability",
    attrib.levels = data.frame(widthvariability_class=
      c("high","moderate","none")),
    u             = c(1,0.4125,0),
    names.u       = c("u.high","u_moderate","u.none"),
    required      = FALSE,
    utility       = FALSE)

# define 1d interpolation end node for bed modification with 
# riprap
# (attribute "bedmodfract_percent" with levels from 0 to 100)

bedmod_riprap <-
  utility.endnode.intpol1d.create(
    name.node   = "bed modification riprap",
    name.attrib = "bedmodfract_percent",
    range       = c(0,100),
    x           = c(0,10,30,100),
    u           = c(1,0.775,0.5625,0.24),
    required    = FALSE,
    utility     = FALSE)

# define 1d interpolation end node for bed modification with 
# other material
# (attribute "bedmodfract_percent" with levels from 0 to 100)

bedmod_other <-
  utility.endnode.intpol1d.create(
    name.node   = "bed modification other",
    name.attrib = "bedmodfract_percent",
    range       = c(0,100),
    x           = c(0,10,30,100),
    u           = c(1,0.775,0.5625,0),
    required    = FALSE,
    utility     = FALSE)

# define combination end node for bed modification
# (attributes "bedmodtype_class" and "bedmodfract_percent")

bedmod <-
  utility.endnode.cond.create(
    name.node     = "bed modification",
    attrib.levels = data.frame(bedmodtype_class=
      c("riprap","other")),
    nodes         = list(bedmod_riprap,bedmod_other),
    required      = FALSE,
    utility       = FALSE)

# define 1d interpolation end node for bank modification with 
# permeable material
# (attribute "bankmodfract_percent" with levels from 0 to 100)

bankmod_perm <-
  utility.endnode.intpol1d.create(
    name.node   = "bank modification perm",
    name.attrib = "bankmodfract_percent",
    range       = c(0,100),
    x           = c(0,10,30,60,100),
    u           = c(1,0.8667,0.675,0.4125,0.24),
    required    = FALSE,
    utility     = FALSE)

# define 1d interpolation end node for bank modification with 
# impermeable material
# (attribute "bankmodfract_percent" with levels from 0 to 100)

bankmod_imperm <-
  utility.endnode.intpol1d.create(
    name.node   = "bank modification imperm",
    name.attrib = "bankmodfract_percent",
    range       = c(0,100),
    x           = c(0,10,30,60,100),
    u           = c(1,0.775,0.5625,0.24,0),
    required    = FALSE,
    utility     = FALSE)

# define combination end node for bank modification
# (attributes "bankmodtype_class" and "bankmodfract_percent")

bankmod <-
  utility.endnode.cond.create(
    name.node     = "bank modification",
    attrib.levels = data.frame(bankmodtype_class=
      c("perm","imperm")),
    nodes         = list(bankmod_perm,bankmod_imperm),
    required      = FALSE,
    utility       = FALSE)

# define 2d interpolation end node for riparian zone width
# (attributes "riparianzonewidth_m" and "riparianzonewidth_m")

riparzone_width <-
  utility.endnode.intpol2d.create(
    name.node   = "riparian zone width",
    name.attrib = c("riverbedwidth_m","riparianzonewidth_m"),
    ranges      = list(c(0,16),c(0,30)),
    isolines    = list(list(x=c(0,16),y=c(0,0)),
                       list(x=c(0,2,10,16),y=c(5,5,15,15)),
                       list(x=c(0,16),y=c(15,15)),
                       list(x=c(0,16),y=c(30,30))),
    u           = c(0.0,0.6,1.0,1.0),
    lead         = 1,
    utility      = FALSE)

# define discrete end node for riparian zone vegetation
# (attriute "riparianzoneveg_class" with levels "natural", 
# "seminatural" and "artificial")

riparzone_veg <-
  utility.endnode.discrete.create(
    name.node     = "riparian zone veg.",
    attrib.levels = data.frame(riparianzoneveg_class=
      c("natural","seminatural","artificial")),
    u             = c(1,0.5625,0),
    required      = FALSE,
    utility       = FALSE)

# define aggregation node for riparian zone

riparzone <-
  utility.aggregation.create(
    name.node = "riparian zone",
    nodes     = list(riparzone_width,riparzone_veg),
    name.fun  = "utility.aggregate.cobbdouglas",
    par       = c(1,1),
    required  = FALSE)

# define aggregation node for ecomorphological state

morphol <-
  utility.aggregation.create(
    name.node = "ecomorphology",
    nodes     = list(widthvar,bedmod,bankmod,riparzone),
    name.fun  = "utility.aggregate.mix",
    par       = c(0.25,0.25,0.25,0.25,0,0,1),
    names.par = c("w_widthvar","w_bedmod","w_bankmod","w_riparzone",
                  "w_add","w_min","w_cobbdouglas"),
    required  = TRUE)

# print individual definitions

print(widthvar)
print(bedmod)

# print all definitions

print(morphol)

# plot objectives hierarchy with attributes

plot(morphol)

# plot individual nodes:

plot(widthvar)
plot(widthvar,par=c(u_moderate=0.2))
plot(bedmod_other)
plot(bankmod)
#plot(riparzone_width)   # too slow for package installation

# plot selected node definitions of a hierarchy

plot(morphol,type="nodes",nodes=c("width variability",
                                  "bed modification other",
                                  "bank modification"))

# evaluate value function for data sets and plot colored hierarchies
# and table

attrib_channelized <- data.frame(widthvariability_class = "none",
                                 bedmodtype_class       = "riprap",
                                 bedmodfract_percent    = 50,
                                 bankmodtype_class      = "imperm",
                                 bankmodfract_percent   = 70,
                                 riverbedwidth_m        = 10,
                                 riparianzonewidth_m    = 5,
                                 riparianzoneveg_class  = "seminatural")
attrib_rehab       <- data.frame(widthvariability_class = "high",
                                 bedmodtype_class       = "riprap",
                                 bedmodfract_percent    = 50,
                                 bankmodtype_class      = "imperm",
                                 bankmodfract_percent   = 20,
                                 riverbedwidth_m        = 15,
                                 riparianzonewidth_m    = 15,
                                 riparianzoneveg_class  = "natural")

res_channelized     <- evaluate(morphol,attrib=attrib_channelized)
res_channelized_add <- evaluate(morphol,attrib=attrib_channelized,
                                par=c(w_add=1,w_min=0,w_cobbdouglas=0))
res_rehab           <- evaluate(morphol,attrib=attrib_rehab)
res_both            <- rbind(res_channelized,res_rehab)
rownames(res_both)  <- c("channelized","rehabilitated")

plot(morphol,u=res_channelized)
plot(morphol,u=res_channelized_add)
plot(morphol,u=res_rehab)
plot(morphol,u=res_rehab,uref=res_channelized)
plot(morphol,u=res_both,type="table",plot.val=FALSE)
plot(morphol,u=res_both,type="table",plot.val=TRUE,print.val=FALSE)
plot(morphol,u=res_both,uref=res_channelized,type="table",plot.val=FALSE)

# consideration of uncertain attribute levels 
# (Higher uncertainty for predicted state after rehabilitation than for 
# observed channelized state.
# Note that the normal distributions lead to a small probability of attribute 
# levels beyond the range for which the value function is defined. This could
# be corrected for by truncating or choosing another distribution. We keep 
# those values to demonstrate that this leads to warnings when evaluating the
# value function for these attribute levels,):

sampsize <- 1000

attrib_channelized_unc <- data.frame(
  widthvariability_class = rep("high",sampsize),
  bedmodtype_class       = rep("riprap",sampsize),
  bedmodfract_percent    = rnorm(sampsize,mean=50,sd=5),
  bankmodtype_class      = rep("imperm",sampsize),
  bankmodfract_percent   = rnorm(sampsize,mean=70,sd=5),
  riverbedwidth_m        = rep(10,sampsize),
  riparianzonewidth_m    = rep(5,sampsize),
  riparianzoneveg_class  = c("seminatural","artificial")[rbinom(sampsize,1,0.5)+1])

attrib_rehab_unc <- data.frame(
  widthvariability_class = c("moderate","high")[rbinom(sampsize,1,0.5)+1],
  bedmodtype_class       = rep("riprap",sampsize),
  bedmodfract_percent    = rnorm(sampsize,mean=50,sd=15),
  bankmodtype_class      = rep("imperm",sampsize),
  bankmodfract_percent   = rnorm(sampsize,mean=20,sd=5),
  riverbedwidth_m        = rnorm(sampsize,mean=10,sd=2),
  riparianzonewidth_m    = rnorm(sampsize,mean=10,sd=2),
  riparianzoneveg_class  = c("natural","seminatural")[rbinom(sampsize,1,0.5)+1])

res_channelized_unc <- evaluate(morphol,attrib=attrib_channelized_unc)
res_rehab_unc       <- evaluate(morphol,attrib=attrib_rehab_unc)

plot(morphol,u=res_channelized_unc)
#plot(morphol,u=res_rehab_unc)
plot(morphol,u=res_rehab_unc,uref=res_channelized_unc)
plot(morphol,u=list(channelized=res_channelized_unc,rehabilitated=res_rehab_unc),
     type="table")
plot(morphol,u=list(channelized=res_channelized_unc,rehabilitated=res_rehab_unc),
     type="table",nodes=c("ecomorphology","riparian zone"))
plot(morphol,u=list(channelized=res_channelized_unc,rehabilitated=res_rehab_unc),
     type="table",levels=2)
plot(morphol,u=list(channelized=res_channelized_unc,rehabilitated=res_rehab_unc),
     uref=res_channelized_unc,
     type="table")
}
