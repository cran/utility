################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.4                                        Peter Reichert 07.01.2017 #
#                                                                              #
################################################################################


# ==============================================================================
# utility node for (potentially) aggregating utility and/or end nodes: 
# class "utility.aggregation"
# ==============================================================================


# constructor:
# ------------

utility.aggregation.create <- 
  function(name.node,          # character(1)
           nodes,              # list of nodes
           name.fun,           # name of aggreg. fun f(u,par)
           par,                # numeric(n)
           names.par    = rep(NA,length(par)),
           required     = FALSE,
           num.required = 1,
           col          = "black",
           shift.levels = 0,
           add.arg.fun  = NULL)
  {
    # consistency checks:
    
    check.ok <- T   
    if ( length(nodes) < 1 )
    {
      cat("*** Warning: No nodes provided","\n")
      check.ok <- F
    }
    utility <- nodes[[1]]$utility
    if ( length(nodes) > 1 )
    {
      for ( i in 2:length(nodes) )
      {
        if ( nodes[[i]]$utility != utility )
        {
          cat("*** Warning: Mixted value and utility nodes",
              "cannot be aggregated","\n")
          check.ok <- F
        }
      }
    }
    if ( ! utility.check.name(name.node,nodes) )
    {
      cat("*** Warning: Node with same name \"",name.node,"\" exists already ",
          "as sub-node","\n")
      check.ok <- F
    }
    if ( ! check.ok )
    {
      cat("*** Warning: Node \"",name.node,"\" could not be constructed","\n",
          sep="")
      return(NA)
    }
    
    # construct class:
    
    node <- list()
    node$name         <- name.node
    node$description  <- "utility/value aggregation node"
    node$type         <- "aggregationnode"
    node$nodes        <- nodes
    node$name.fun     <- name.fun
    node$par          <- par
    node$names.par    <- names.par
    node$required     <- required
    node$num.required <- num.required
    node$utility      <- utility
    node$col          <- col
    node$shift.levels <- shift.levels
    node$add.arg.fun  <- add.arg.fun
    class(node)       <- "utility.aggregation" 
    
    # return class
    
    #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
    return(node)
  }


# update parameter values:
# ------------------------

updatepar.utility.aggregation <- function(x,par=NA,...)
{
  node <- x
  
  # check availabiliy of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # update conditional nodes:
  
  n <- node
  for ( i in 1:length(n$par) )
  {
    if ( ! is.na(n$names.par[i]) )
    {
      ind <- which(n$names.par[i] == names(par) )
      if ( length(ind) > 1 )
      {
        warning("Node \"",node$name,"\": multiple occurrences of parameter",
                names(par)[ind[1]],sep="")
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$par[i] <- par[ind]
      }
    } 
  }
  for ( i in 1:length(n$nodes) )
  {
    n$nodes[[i]] <- updatepar(n$nodes[[i]],par)
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.aggregation <- function(x,
                                         attrib,   # data.frame
                                         par=NA,
                                         ...)
{
  node <- x
  
  # check input:
  
  if ( ! is.data.frame(attrib) )
  {
    warning("Node \"",node$name,"\": attrib must be a data frame",sep="")
    return(NA)
  }
  
  # update parameters:
  
  n <- updatepar(node,par)
  
  # evaluate nodes:
  
  u <- evaluate(n$nodes[[1]],attrib)
  ind <- !is.na(u) & (u<0 | u>1)
  if ( sum(ind) > 0 )
  {
    warning("Node \"",node$name,"\": node \"",n$nodes[[1]]$name,"\" produced values outside [0,1]: ",
            paste(u[ind],collapse=","),sep="")
  }
  if ( ! is.data.frame(u) )
  {
    u <- as.data.frame(u)
    names(u) <- n$nodes[[1]]$name
  }
  required <- n$nodes[[1]]$required
  nodenames <- n$nodes[[1]]$name
  if ( length(n$nodes) > 1 )
  {
    for ( i in 2:length(n$nodes) )
    {
      u.i <- evaluate(n$nodes[[i]],attrib)
      ind <- !is.na(u) & (u<0 | u>1)
      if ( sum(ind) > 0 )
      {
        warning("Node \"",node$name,"\": node \"",n$nodes[[i]]$name,"\" produced values outside [0,1]: ",
                paste(u.i[ind],collapse=","),sep="")
      }
      if ( ! is.data.frame(u.i) )
      {
        u.i <- as.data.frame(u.i)
        names(u.i) <- n$nodes[[i]]$name
      }
      u <- cbind(u,u.i)
      nodenames[i] <- n$nodes[[i]]$name
      required[i]  <- n$nodes[[i]]$required 
    }
  }
  if ( length(unique(nodenames)) != length(nodenames) )
  {
    warning("Node \"",node$name,"\": node names are not unique:",
            paste(nodenames,collapse=","))
    u.agg <- as.data.frame(rep(NA,nrow(attrib)))
    names(u.agg) <- n$name
    u <- cbind(u.agg,u)
    rownames(u) <- rownames(attrib)
    return(u)
  }   
  
  # return results:
  
  u.agg.input <- as.matrix(u[,nodenames])
  if ( length(n$add.arg.fun) > 0 )
  {
    u.agg <- apply(u.agg.input,1,n$name.fun,n$par,n$add.arg.fun)
  }
  else
  {
    u.agg <- apply(u.agg.input,1,n$name.fun,n$par)
  }
  res.ok <- apply(u.agg.input,1,utility.check.required,
                  required,n$num.required)
  u.agg <- ifelse(res.ok,u.agg,NA)
  u.agg <- as.data.frame(u.agg)
  names(u.agg) <- n$name
  ind <- !is.na(u.agg) & (u.agg<0 | u.agg>1)
  if ( sum(ind)  > 0 )
  {
    warning("Node \"",node$name,"\": aggregation technique \"",n$name.fun,"\" produced values outside of [0,1]: ",
            paste(u.agg[ind],collapse=","),sep="")
  }
  u <- cbind(u.agg,u)
  rownames(u) <- rownames(attrib)
  
  return(u)
}


# print:
# -----

print.utility.aggregation <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.aggregation <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  for ( i in 1:length(node$nodes) )
  {
    string1 <- "nodes:          "
    if ( i > 1 ) string1 <- "                "
    string2 <- node$nodes[[i]]$name
    if ( node$nodes[[i]]$type == "endnode" ) 
    {
      num.space <- max(1,15-nchar(node$nodes[[i]]$name))
      string2 <- paste(string2,
                       paste(rep(" ",num.space),collapse=""),
                       "(end node)",sep="") 
    }     
    cat(string1,string2,"\n")
  }
  cat("function:       ",node$name.fun,"\n")
  names.par <- ifelse(is.na(node$names.par),"",node$names.par)
  cat("parameters:","\n")
  print(data.frame(names.par=names.par,par=node$par))
  if ( length(node$add.arg.fun) > 0 ) print(node$add.arg.fun)
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("required nodes: ",node$num.required,"\n")
  for ( i in 1:length(node$nodes) ) 
  {
    cat("***","\n")
    summary(node$nodes[[i]])
  }
}


# plot:
# -----

plot.utility.aggregation <- 
  function(x,
           u           = NA,
           uref        = NA,
           par         = NA,
           type        = c("hierarchy","table","node","nodes"),
           nodes       = NA,
           col         = utility.calc.colors(),
           gridlines   = c(0.2,0.4,0.6,0.8),
           main        = "",
           cex.main    = 1,
           cex.nodes   = 1,
           cex.attrib  = 1,
           f.reaches   = 0.2,
           f.nodes     = 0.2,
           with.attrib = TRUE,
           levels      = NA,
           plot.val    = TRUE,
           print.val   = TRUE,
           two.lines   = FALSE,
           ...)
{
    node <- x
    n <- updatepar(node,par)
    utility.plot(node        = n,
                 u           = u,
                 uref        = uref,
                 type        = type,
                 nodes       = nodes,
                 col         = col,
                 gridlines   = gridlines,
                 main        = main,
                 cex.main    = cex.main,
                 cex.nodes   = cex.nodes,
                 cex.attrib  = cex.attrib,
                 f.reaches   = f.reaches,
                 f.nodes     = f.nodes,
                 with.attrib = with.attrib,
                 levels      = levels,
                 plot.val    = plot.val,
                 print.val   = print.val,
                 two.lines   = two.lines,
                 ...)
  }


