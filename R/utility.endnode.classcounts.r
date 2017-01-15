################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.4                                        Peter Reichert 05.09.2016 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for valuing counts in discrete classes: 
# class "utility.endnode.classcounts"
# ==============================================================================


# constructor:
# ------------

utility.endnode.classcounts.create <- function(name.node,        # character(1)
                                               name.attrib,      # character(n)
                                               u.max.inc,        # list (n) of vect (>=1)
                                               names.u.max.inc = list(),
                                               exceed.next     = TRUE,
                                               utility         = TRUE,
                                               required        = FALSE,
                                               col             = "black",
                                               shift.levels    = 0)
{
  # consistency checks:
  
  check.ok <- T   
  n <- length(name.attrib)
  if ( length(u.max.inc) != n )
  {
    cat("*** Warning: Number of elements of u.max.inc not equal to number of elements of name.attrib:",
        length(u.max.inc),n,"\n")
    check.ok <- F
  }
  for ( i in 1:n )
  {
    if ( !is.vector(u.max.inc[[i]]) )
    {
      cat("*** Warning: Eelements of u.max.inc must be vectors","\n")
      check.ok <- F
    }
  }
  if ( length(names.u.max.inc) != 0 & length(names.u.max.inc) != n )
  {
    cat("*** Warning: Number of elements of names.u.max.inc not equal to zero or to the number of elements of name.attrib:",
        length(names.u.max.inc),n,"\n")
    check.ok <- F
  }
  if ( ! check.ok )
  {
    cat("*** Warning: node \"",name.node,"\" could not be constructed","\n",
        sep="")
    return(NA)
  }
  
  # construct class:
  
  node <- list()
  node$name            <- name.node
  node$description     <- "utility/value class counts end node"
  node$type            <- "endnode"
  node$attrib          <- name.attrib
  node$u.max.inc       <- u.max.inc
  for ( i in 1:n )
  {
    l <- length(node$u.max.inc[[i]])
    if ( l < n+2-i ) node$u.max.inc[[i]] <- c(node$u.max.inc[[i]],rep(0,n+2-i-l))
    if ( l > n+2-i ) node$u.max.inc[[i]] <- node$u.max.inc[[i]][1:(n+2-i)]
  }
  if ( length(node$names.u.max.inc) ==  n )
  {
    for ( i in 1:n )
    {
      l <- length(node$names.u.max.inc[[i]])
      if ( l < n+2-i ) node$names.u.max.inc[[i]] <- c(node$names.u.max.inc[[i]],rep(NA,n+2-i-l))
      if ( l > n+2-i ) node$names.u.max.inc[[i]] <- node$names.u.max.inc[[i]][1:(n+2-i)]
    }
  }
  node$names.u.max.inc <- names.u.max.inc
  node$exceed.next     <- exceed.next
  node$required        <- required
  node$utility         <- utility
  node$col             <- col
  node$shift.levels    <- shift.levels
  class(node)          <- "utility.endnode.classcounts" 
  
  # print and return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.classcounts <- function(x,par=NA,...)
{
  node <- x
  
  n <- length(node$attrib)
  
  # check availability of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # check availability of parameter names
  
  if ( length(node$names.u.max.inc) !=  n ) return(node)
  
  # update adequate values:
  
  for ( i in 1:length(node$attrib) )
  {
    for ( j in 1:(n+2-i) )
    {
      if ( ! is.na(node$names.u.max.inc[[i]][j]) )
      {
        ind <- which(node$names.u.max.inc[[i]][j] == names(par) )
        if ( length(ind) > 1 )
        {
          warning("Node \"",node$name,"\": multiple occurrences of parameter",
                  names(par)[ind[1]])
          ind <- ind[1]
        }
        if ( length(ind) == 1 )
        {
          node$u.max.inc[[i]][j] <- par[ind]
        }
      }
    }
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.classcounts <- function(x,
                                                 attrib,   # data.frame
                                                 par = NA,
                                                 ...)
{
  node <- x
  n <- length(node$attrib)
  
  # update parameters:
  
  node <- updatepar(node,par)

  # extract attributes:
  
  if ( is.data.frame(attrib) | is.matrix(attrib) )
  {
    ind <- match(node$attrib,colnames(attrib))
    if ( sum(ifelse(is.na(ind),1,0)) > 0 )
    {
      warning("Node \"",node$name,"\": attribute(s) \"",
              paste(node$attrib[is.na(ind)],collapse=","),"\" not found",sep="")
      return(rep(NA,nrow(attrib)))
    }
    a <- attrib[,ind]
  }
  else
  {
    if ( ! is.vector(attrib) )
    {
      warning("Node \"",node$name,"\": unknown format of attribute(s) \"",node$attrib,"\"",sep="")
      return(NA)
    }
    if ( length(names(attrib)) == 0 )
    {
      if ( length(attrib) == 2 )
        a <- as.matrix(attrib,nrow=1)
    }
    else
    {
      ind <- match(node$attrib,names(attrib))
      if ( sum(ifelse(is.na(ind),1,0)) > 0 )
      {
        warning("Node \"",node$name,"\": attribute(s) \"",
                paste(node$attrib[is.na(ind)],collapse=","),"\" not found",sep="")
        return(rep(NA,nrow(attrib)))
      }
      a <- as.matrix(attrib[ind],nrow=1)
    }
  }
  
  # evaluate results:
  
  u <- rep(NA,nrow(a))
  for ( k in 1:nrow(a) )
  {
    att <- as.numeric(a[k,])
    i <- match(TRUE,att>0)
    if ( is.na(i) )
    {
      if ( sum(!is.na(att)) > 0 ) u[k] <- 0
    }
    else
    {
      # basic value:
      u[k] <- node$u.max.inc[[i]][1]
      
      # increment for multiplicities at the maximum level:
      u[k] <- u[k] + (att[i]-1)*node$u.max.inc[[i]][2]
      
      # increment for multiplicities at lower levels:
      if ( i < n )
      {
        for ( j in 1:(n-i) ) u[k] <- u[k] + att[i+j]*node$u.max.inc[[i]][2+j]
      }
      
      # check maximum:
      u[k] <- min(1,u[k])
      if ( i > 1 & !node$exceed.next )
      {
        u[k] <- min(node$u.max.inc[[i-1]][1],u[k])
      }
    }
  }
  
  # return results:
  
  return(u)
}


# print:
# -----

print.utility.endnode.classcounts <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.classcounts <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  cat("attribute(s):   ",paste(node$attrib,collapse=","),"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("basic level, multiplicity increments","\n")
  for ( i in 1:length(node$attrib) )
  {
    cat(paste(node$u.max.inc[[i]],collapse=", "),"\n")
  }
}


# plot:
# -----

plot.utility.endnode.classcounts <- 
  function(x,
           par       = NA,
           col       = utility.calc.colors(),
           gridlines = c(0.2,0.4,0.6,0.8),
           main      = "",
           cex.main  = 1,
           ...)
  {
    # plot frame:
    
    node <- x
    space <- 0.2
    n <- updatepar(node,par)
    title <- main; if ( nchar(title) == 0 ) title <- n$name
    funtype <- "utility"; if ( !n$utility ) funtype <- "value"
    n.attrib <- length(n$attrib)
    u.max.inc <- matrix(0,nrow=n.attrib+1,ncol=n.attrib)
    colnames(u.max.inc) <- n$attrib
    for ( i in 1:n.attrib ) u.max.inc[1:(n.attrib+2-i),i] <- node$u.max.inc[[i]]
    print(u.max.inc)
    barplot(u.max.inc,main=title,ylab=paste(funtype,"(base + inc)"),
            cex.main=cex.main,xlim=c(0,n.attrib*(1+space)),ylim=c(0,1),
            space=space,beside=FALSE,xaxs="i",yaxs="i")
    max.val <- 0.995*rep(1,n.attrib)
    if ( !n$exceed.next ) max.val <- c(0.995,u.max.inc[1,1:(n.attrib-1)])
    for ( i in 1:n.attrib )
    {
      lines(space+(i-1)*(1+space)+0.5+c(-0.5,0.5),max.val[i]*c(1,1),col="red",lwd=2)
    }
  }



# test code:

# library(utility)
# 
# n <- utility.endnode.classcounts.create(
#   name.node = "test",
#   name.attrib = c("a","b","c"),
#   u.max       = list(c(0.8,0.05,0.01),
#                      c(0.5,0.05,0.01),
#                      c(0.0,0.01)),
#   exceed.next = FALSE,
#   utility     = FALSE)
# 
# attrib <- data.frame(a=c(0,0,0,1,2,3,0),
#                      b=c(1,2,9,2,2,2,0),
#                      c=c(5,9,2,8,7,8,0))
# 
# values <- evaluate(n,attrib)
# 
# print(n)
# plot(n)
# print(attrib)
# print(values)

