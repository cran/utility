\name{utility.aggregation.create}

\alias{utility.aggregation.create}

\title{Construct an aggregation node}

\description{Function to construct an aggregation node for value or utlity functions.}

\usage{
utility.aggregation.create(name.node, 
                           nodes, 
                           name.fun, 
                           par, 
                           names.par    = rep(NA, length(par)), 
                           required     = FALSE, 
                           num.required = 1, 
                           col          = "black", 
                           shift.levels = 0,
                           add.arg.fun  = NULL)
}

\arguments{
  \item{name.node}{
name of the node to be constructed as a character string.
}
  \item{nodes}{
list of nodes to be aggregated.
}
  \item{name.fun}{
name of the function to be used for aggregation.
This function must accept the arguments \code{u} and \code{par} which
pass a vector of values or utilities to be aggregated and the parameters
of the function, respectively.
The function can have an additional argument specified below as \code{add.arg.fun}.
The function must then return the corresponding aggregated value or utility.\cr
Examples of functions provided by the package are \cr
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
Follow the links for the aggregation functions for their use, for the underlying mathematical expressions and for graphical illustrations.
}                                                          
  \item{par}{
numeric vector of parameter values to be passed to the function specified unter \code{name.fun}.
}
  \item{names.par}{
(optional) vector of parameter names corresponging to the vector of values specified under \code{par}.
Only required to provide access to the values through a named parameter vector.
}
  \item{required}{
(optional) logical variable indicating if the value of this node is required for aggregation at the next higher level.
If this variable is \code{TRUE}, aggregation at the next higher level is not possible if this node returns NA.
Default value is \code{FALSE}.
}
  \item{num.required}{
number of lower-level values or utilities that must at least be available to make the evaluation possible.
}
  \item{col}{
(optional) color used for plotting the bounding box of the node in the objective hierarchy.
Default value is \code{"black"}.
}
  \item{shift.levels}{
(optional) number of hierarchical levels by which the node in the objective hierarchy is shifted to make a branch fit better to other branches.
Default value is \code{0}.
}
  \item{add.arg.fun}{
(optional) an additional argument to the aggregation function \code{name.fun}.
The value(s) given here will alway be passed to the aggregation function.
}
}

\value{
The function returns the created object of type \code{utility.aggregation} with the properties specified in the arguments of the function.
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
Print, evaluate and plot the node with \cr\cr
\code{\link{print.utility.aggregation}}, \cr
\code{\link{summary.utility.aggregation}}, \cr
\code{\link{evaluate.utility.aggregation}} and \cr
\code{\link{plot.utility.aggregation}}. \cr\cr
Create end nodes with \cr\cr
\code{\link{utility.endnode.discrete.create}}, \cr
\code{\link{utility.endnode.intpol1d.create}}, \cr
\code{\link{utility.endnode.intpol2d.create}}, \cr
\code{\link{utility.endnode.parfun1d.create}}, \cr
\code{\link{utility.endnode.cond.create}}, or \cr
\code{\link{utility.endnode.firstavail.create}}. \cr\cr
Create conversion nodes with \cr\cr
\code{\link{utility.conversion.intpol.create}}, or \cr
\code{\link{utility.conversion.parfun.create}}.
}

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
  #plot(riparzone_width)
  
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
  plot(morphol,u=res_both,type="table")

  # consideration of uncertain attribute levels (higher uncertainty for 
  # predicted state after rehabilitation than for observed channelized state):
  
  sampsize <- 2000
  
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
  plot(morphol,u=res_rehab_unc)
  plot(morphol,u=res_rehab_unc,uref=res_channelized_unc)
}

