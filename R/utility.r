################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.01                                       Peter Reichert 18.10.2012 #
#                                                                              #
################################################################################


# ==============================================================================
# registration of member functions
# ==============================================================================


updatepar     <- function(x, ...) UseMethod("updatepar")
evaluate      <- function(x, ...) UseMethod("evaluate")
evaluate.cond <- function(x, ...) UseMethod("evaluate.cond")

# in addition, we support the functions plot, print and summary


# ==============================================================================
# simple parametric utility functions
# ==============================================================================


utility.fun.exp <- function(attrib,par)   # par[1]:  absolute risk aversion
{                                         # par[2]:  minimum of attribute range (default=0)
                                          # par[3]:  maximum of attribute range (default=1)
   atrans <- attrib
   if ( length(par) >= 3 ) atrans <- (attrib-par[2])/(par[3]-par[2])
   if ( par[1] == 0 ) return(atrans)
   return((1-exp(-atrans*par[1]))/(1-exp(-par[1])))
}


# ==============================================================================
# utility aggregation functions
# ==============================================================================


utility.aggregate.add <- function(u,par)  # par[i]: weight of u[i]
{
   # check input:
   
   if ( length(u) != length(par) )
   {
      warning("Length of utilities/values and weights not equal:",
              length(u),length(par))
      return(NA)
   }
   ind <- which(!is.na(u))
   if ( length(ind) == 0 ) return(NA)
   if ( sum( par < 0 ) > 0 )
   {
      warning("Parameter of additive aggregation smaller than zero")
      return(NA)
   }
   
   # calculate aggregated value

   s <- sum(par[ind])
   if ( s <= 0 ) return(NA)
   u.agg <- sum(par[ind]*u[ind])/s
   
   return(as.numeric(u.agg))
}


utility.aggregate.min <- function(u,par=NA)
{
   # check input:
   
   ind <- which(!is.na(u))
   if ( length(ind) == 0 ) return(NA)
   
   # calculate aggregated value
   
   u.agg <- min(u[ind])
   
   return(as.numeric(u.agg))
}


utility.aggregate.max <- function(u,par=NA)
{
   # check input:
   
   ind <- which(!is.na(u))
   if ( length(ind) == 0 ) return(NA)
   
   # calculate aggregated value
   
   u.agg <- max(u[ind])
   
   return(as.numeric(u.agg))
}


utility.aggregate.mult <- function(u,par) 
{
   # check input:
   
   if ( length(u) != length(par) )
   {
      warning("Length of utilities/values and weights not equal:",
              length(u),length(par))
      return(NA)
   }
   ind <- which(!is.na(u))
   if ( length(ind) == 0 ) return(NA)
   if ( length(ind) == 1 )
   {
      return(as.numeric(u[ind]))
   }
   if ( sum( par < 0 | par > 1 ) > 0 )
   {
      warning("Parameter of multiplicative aggregation",
              "smaller than zero or larger than unity")
      return(NA)
   }
   
   # define numerical parameter:
   
   eps <- 1e-3   # maximum deviation of sum(par) from unity to use additive fcn 

   # rescale weights:

   s <- sum(par)   
   fact <- s/sum(par[ind])
   ki <- fact*par[ind]
      
   # calculate additive utility function if sum close to unity:

   if ( s > 1-eps & s < 1+eps )
   {
      return(utility.aggregate.add(u,par))
   }
   
   # calculate multiplicative utility function if sum not close to unity:
   
   # calculate k: 
   # (Keeney and Raiffa, Decisions with multiple objectives, 1976,
   # pp. 307, 347-348)
   
   if ( s < 1 )
   {
      lower <- 1
      i <- 0
      while ( utility.aggregate.mult.f.root(lower,ki) < 0 )
      {
         lower <- 0.1*lower
         i <- i+1
         if ( i > 20 )
         {
            warning("Problem solving equation for scaling constant")
            return(NA)
         }
      }
      upper <- 1
      i <- 0
      while ( utility.aggregate.mult.f.root(upper,ki) > 0 )
      {
         upper <- 10*upper
         i <- i+1
         if ( i > 20 )
         {
            warning("Problem solving equation for scaling constant")
            return(NA)
         }
      }
      k <- uniroot(utility.aggregate.mult.f.root,ki=ki,
                   lower=lower,upper=upper)$root
   }
   else  # s > 1
   {
      upper <- -0.1
      i <- 0
      while ( utility.aggregate.mult.f.root(upper,ki) < 0 )
      {
         upper <- 0.1*upper
         i <- i+1
         if ( i > 20 )
         {
            warning("Problem solving equation for scaling constant")
            return(NA)
         }
      }
      k <- uniroot(utility.aggregate.mult.f.root,ki=ki,
                   lower=-1,upper=upper)$root 
   }

   # evaluate multiplicative utility function:
   
   u.agg <- 1  
   for ( i in 1:length(ki) )
   {
      if ( !is.na(u[ind][i]) ) u.agg <- u.agg * (k*ki[i]*u[ind][i]+1) 
   }
   u.agg <- (u.agg - 1)/k
      
   return(as.numeric(u.agg))
}


utility.aggregate.mult.f.root <- function(k,ki)
{
   res <- 1
   for ( i in 1:length(ki) )
   {
      res <- res * ( 1 + k * ki[i] )
   }
   res <- 1 + k - res
   return(res)
}


utility.aggregate.cobbdouglas <- function(u,par) 
{
   # check input:
   
   if ( length(u) != length(par) )
   {
      warning("Length of utilities/values and weights not equal:",
              length(u),length(par))
      return(NA)
   }
   ind <- which(!is.na(u))
   if ( length(ind) == 0 ) return(NA)
   if ( sum( par < 0 ) > 0 )
   {
      warning("Parameter of Cobb-Douglas aggregation smaller than zero")
      return(NA)
   }
   
   # calculate aggregated value

   s <- sum(par[ind])
   if ( s <= 0 ) return(NA)
   u.agg <- 1
   for ( i in 1:length(ind) )
   {
      if ( par[ind][i]>0 ) u.agg <- u.agg*u[ind][i]^(par[ind][i]/s)
   }
   
   return(as.numeric(u.agg))
}


utility.aggregate.mix <- function(u,par)  # par[i]: weight of u[i]
{                                         # par[n+j]: weight of technique j
   # check input:                         # (j = add, min, cobbdouglas)
   
   n <- length(u)
   if ( n+3 != length(par) )
   {
      warning("Length of parameter vector must be equal to",
             "length of utilities/values plus three:",
             length(par),length(u))
      return(NA)
   }
   s <- sum(par[n+(1:3)])
   if ( s <= 0 | sum(par[n+(1:3)]<0) > 0 )
   {
      warning("Weights of aggregation techniques to average",
              "cannot be negative or not all of them equal to zero")
      return(NA)
   }
   
   u.add         <- utility.aggregate.add(u,par[1:n])
   u.min         <- utility.aggregate.min(u)
   u.cobbdouglas <- utility.aggregate.cobbdouglas(u,par[1:n])
   
   if ( is.na(u.add) | is.na(u.min) | is.na(u.cobbdouglas) ) return(NA)
   u.agg <- (par[n+1]*u.add + par[n+2]*u.min + par[n+3]*u.cobbdouglas)/s

   return(u.agg)
}


# ==============================================================================
# auxiliary functions
# ==============================================================================


utility.calc.colors <- function(n=5)
{
   if ( n < 2 ) return("black")
   if ( n < 3 ) return(c("tomato","blue"))
   if ( n < 4 ) return(c("tomato","yellow","blue"))
   if ( n < 5 ) return(c("tomato","yellow","green","blue"))
   if ( n < 6 ) return(c("tomato","orange","yellow","lightgreen","lightblue"))

   red    <- col2rgb("tomato")/255
   orange <- col2rgb("orange")/255
   yellow <- col2rgb("yellow")/255
   green  <- col2rgb("lightgreen")/255
   blue   <- col2rgb("lightblue")/255
   red.orange    <- (2*red+orange)/3
   orange.red    <- (red+2*orange)/3
   orange.yellow <- (2*orange+yellow)/3
   yellow.orange <- (orange+2*yellow)/3
   yellow.green  <- (2*yellow+green)/3
   green.yellow  <- (yellow+2*green)/3
   green.blue    <- (2*green+blue)/3
   blue.green    <- (1.5*green+blue)/2.5

   u <- (1:n)/(n+1)
   cols <- rep(NA,n)

   for ( i in 1:length(u) )
   {
      if( u[i]<0.2 )
      { 
         col <- (1-u[i]/0.2) * red+
                u[i]/0.2 * red.orange
      }
      if( 0.2<=u[i] & u[i]<0.4 )
      {
         col <- (1-(u[i]-0.2)/0.2) * orange.red +
                (u[i]-0.2)/0.2 * orange.yellow 
      }
      if( 0.4<=u[i] & u[i]<0.6 ) 
      {
         col <- (1-(u[i]-0.4)/0.2) * yellow.orange +
                (u[i]-0.4)/0.2 * yellow.green
      }
      if( 0.6<=u[i] & u[i]<0.8 ) 
      {
         col <- (1-(u[i]-0.6)/0.2) * green.yellow +
                (u[i]-0.6)/0.2 * green.blue
      }
      if( 0.8<=u[i] ) 
      {
         col <- (1-(u[i]-0.8)/0.2) * blue.green +
                (u[i]-0.8)/0.2 * blue
      }
      cols[i] <- rgb(col[1],col[2],col[3])
   }
   return(cols)
}


utility.get.colors <- function(u,col=utility.calc.colors())
{
   col.ind <- 1 + floor(u*length(col)*0.99999)
   cols <- col[col.ind]
   cols <- ifelse(is.na(col.ind),"white",cols)
   return(cols)
}


utility.get_y_belowandabove <- function(x,y,xout,yref)
{
   y.res <- c(below=NA,above=NA) 
   if ( xout<min(x) | xout>max(x) ) return(y.res)
   x.lower <- x[-length(x)]
   x.upper <- x[-1]
   ind <- which(ifelse( (xout>=x.lower & xout<=x.upper) |
                        (xout<=x.lower & xout>=x.upper) ,T,F ))
   if ( length(ind) == 0 ) return(y.res)
   y.vals <- rep(NA,length(ind))
   for ( i in 1:length(ind) )
   {
      if ( x[ind[i]+1] == x[ind[i]] )
      {
         if ( (y[ind[i]]>yref) & (y[ind[i]+1]>yref) )
         {
            y.vals[i] <- min(y[ind[i]],y[ind[i]+1])
         }
         else
         {
            if ( (y[ind[i]]<yref) & (y[ind[i]+1]<yref) )
            {
               y.vals[i] <- max(y[ind[i]],y[ind[i]+1])
            }
            else
            {
               y.vals[i] <- yref
            }
         } 
      }
      else
      {
         y.vals[i] <- y[ind[i]] + (xout-x[ind[i]])/(x[ind[i]+1]-x[ind[i]])*
                                  (y[ind[i]+1]-y[ind[i]])
      }
   }
   if ( sum(y.vals<=yref) > 0 ) y.res["below"] <- max(y.vals[y.vals<=yref])
   if ( sum(y.vals>=yref) > 0 ) y.res["above"] <- min(y.vals[y.vals>=yref])
   return(y.res)
}


utility.intpol.multiple <- function(x,xs,ys)
{
   ind <- !is.na(xs) & !is.na(ys)                                 
   if ( sum(ind) < 2 ) return(NA)
   xs.loc <- xs[ind]
   ys.loc <- ys[ind]
   
   ind.below <- which(xs.loc<=x)
   if ( length(ind.below) == 0 ) return(NA)
   ind.above <- which(xs.loc>=x)
   if ( length(ind.above) == 0 ) return(NA)  
   xs.below <- xs.loc[ind.below]
   ys.below <- ys.loc[ind.below]
   xs.above <- xs.loc[ind.above]
   ys.above <- ys.loc[ind.above]
   
   ind.max.below <- which.max(xs.below)
   x.below <- xs.below[ind.max.below]
   y.below <- ys.below[ind.max.below]
   ind.min.above <- which.min(xs.above)
   x.above <- xs.above[ind.min.above]
   y.above <- ys.above[ind.min.above]

   if ( x.above == x.below )
   {
      y <- mean(y.above,y.below)
   }
   else
   {   
      y <- ( y.above*(x-x.below) + y.below*(x.above-x) ) / (x.above-x.below)
   }
   
   return(y)
}


utility.intpol2d <- function(xy,isolines,levels,lead=0)
{
   ind <- order(levels)
   z <- apply(xy,1,utility.intpol2d.pair,isolines[ind],levels[ind],lead)

   return(z)
}


utility.intpol2d.pair <- function(xy,isolines,levels,lead=0)
{
   # initialize u:
   
   z <- rep(NA,2)
   nam <- c("x","y")
   
   xy <- as.numeric(xy)
   if( is.na(xy[1]) | is.na(xy[2]) ) return(NA)
   
   for ( lead.current in 1:2 )
   {
      ind.x <- lead.current
      ind.y <- 3-ind.x
      nam.x <- nam[ind.x]
      nam.y <- nam[ind.y]
      if ( lead == 0 | lead == ind.x )
      {
         for ( i in 2:length(isolines) )
         {
            n.1 <- length(isolines[[i-1]][[nam.x]])
            n.2 <- length(isolines[[i]][[nam.x]])
            if ( xy[ind.x] >= min(isolines[[i-1]][[nam.x]]) & 
                 xy[ind.x] <= max(isolines[[i-1]][[nam.x]]) )
            {
               y.1 <- utility.get_y_belowandabove(x = isolines[[i-1]][[nam.x]],
                                                  y = isolines[[i-1]][[nam.y]],
                                                  xout = xy[ind.x],
                                                  yref = xy[ind.y])
               if ( xy[ind.x] >= min(isolines[[i]][[nam.x]]) & 
                    xy[ind.x] <= max(isolines[[i]][[nam.x]]) )
               {
                  # x coordinate of xy intersects contour lines at
                  # levels i-1 and i
            
                  y.2 <- utility.get_y_belowandabove(x = isolines[[i]][[nam.x]],
                                                     y = isolines[[i]][[nam.y]],
                                                     xout = xy[ind.x],
                                                     yref = xy[ind.y])
                  val <- utility.intpol.multiple(x  = xy[ind.y],
                                                 xs = c(y.1,y.2),
                                                 ys = c(rep(levels[i-1],2),
                                                              rep(levels[i],2)))
                  if ( ! is.na(val) )
                  {
                     z[lead.current] <- val
                     break
                  }
               }
               else  # within range of line at level i-1, 
                     # outside of range at level i
               {
                  if ( xy[ind.x] > max(isolines[[i]][[nam.x]]) )
                  { 
                     # x coordinate of xy intersects contour line at
                     # level i-1 but is larger than maximum x at level i
                     
                     ratio.1 <- NA
                     y.2.1 <- NA
                     z.2.1 <- NA
                     if ( xy[ind.x] < isolines[[i-1]][[nam.x]][1] )
                     {
                        ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                                   (isolines[[i]][[nam.x]][1]-
                                                isolines[[i-1]][[nam.x]][1])
                        y.2.1 <- isolines[[i-1]][[nam.y]][1] +
                                 ratio.1*(isolines[[i]][[nam.y]][1]-
                                                isolines[[i-1]][[nam.y]][1])
                        z.2.1 <- levels[[i-1]] +
                                 ratio.1*(levels[[i]]-levels[[i-1]])
                     }
                     ratio.n <- NA
                     y.2.n <- NA
                     z.2.n <- NA
                     if ( xy[ind.x] < isolines[[i-1]][[nam.x]][n.1] )
                     {
                        ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                                   (isolines[[i-1]][[nam.x]][n.1]-
                                                isolines[[i]][[nam.x]][n.2])
                        y.2.n <- isolines[[i-1]][[nam.y]][n.1] +
                                 ratio.n*(isolines[[i]][[nam.y]][n.2]-
                                                isolines[[i-1]][[nam.y]][n.1])
                        z.2.n <- levels[[i-1]] +
                                 ratio.n*(levels[[i]]-levels[[i-1]])
                     }
                     val <- utility.intpol.multiple(x  = xy[ind.y],
                                                    xs = c(y.1,y.2.1,y.2.n),
                                                    ys = c(rep(levels[i-1],2),
                                                                  z.2.1,z.2.n))
                     if ( ! is.na(val) )
                     {
                        z[lead.current] <- val
                        break
                     }
                  }
                  else # xy[ind.x] < min(isolines[[i]][[nam.x]])
                  {
                     # x coordinate of xy intersects contour line
                     # at level i-1 but is smaller than minimum x at level i

                     ratio.1 <- NA
                     y.2.1 <- NA
                     z.2.1 <- NA
                     if ( xy[ind.x] > isolines[[i-1]][[nam.x]][1] )
                     {
                        ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                                   (isolines[[i]][[nam.x]][1]-
                                                isolines[[i-1]][[nam.x]][1])
                        y.2.1 <- isolines[[i-1]][[nam.y]][1] +
                                 ratio.1*(isolines[[i]][[nam.y]][1]-
                                                isolines[[i-1]][[nam.y]][1])
                        z.2.1 <- levels[[i-1]] +
                                 ratio.1*(levels[[i]]-levels[[i-1]])
                     }
                     ratio.n <- NA
                     y.2.n <- NA
                     z.2.n <- NA
                     if ( xy[ind.x] > isolines[[i-1]][[nam.x]][n.1] )
                     {
                        ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                                   (isolines[[i-1]][[nam.x]][n.1]-
                                                isolines[[i]][[nam.x]][n.2])
                        y.2.n <- isolines[[i-1]][[nam.y]][n.1] +
                                 ratio.n*(isolines[[i]][[nam.y]][n.2]-
                                                isolines[[i-1]][[nam.y]][n.1])
                        z.2.n <- levels[[i-1]] +
                                 ratio.n*(levels[[i]]-levels[[i-1]])
                     }
                     val <- utility.intpol.multiple(x  = xy[ind.y],
                                                    xs = c(y.1,y.2.1,y.2.n),
                                                    ys = c(rep(levels[i-1],2),
                                                                  z.2.1,z.2.n))
                     if ( ! is.na(val) )
                     {
                        z[lead.current] <- val
                        break
                     }
                  }
               }
            }
            else  # outside of range of line at level i-1
            {
               if ( xy[ind.x] >= min(isolines[[i]][[nam.x]]) & 
                    xy[ind.x] <= max(isolines[[i]][[nam.x]]) )
               {
                  y.2 <- utility.get_y_belowandabove(x = isolines[[i]][[nam.x]],
                                                     y = isolines[[i]][[nam.y]],
                                                     xout = xy[ind.x],
                                                     yref = xy[ind.y])                  

                  if ( xy[ind.x] > max(isolines[[i-1]][[nam.x]]) )
                  { 
                     # x coordinate of xy intersects isoline 
                     # at level i but is larger than maximum x at level i-1

                     ratio.1 <- NA
                     y.1.1 <- NA
                     z.1.1 <- NA
                     if ( xy[ind.x] < isolines[[i]][[nam.x]][1] )
                     {
                        ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                                   (isolines[[i]][[nam.x]][1]-
                                                isolines[[i-1]][[nam.x]][1])
                        y.1.1 <- isolines[[i-1]][[nam.y]][1] +
                                 ratio.1*(isolines[[i]][[nam.y]][1]-
                                                isolines[[i-1]][[nam.y]][1])
                        z.1.1 <- levels[[i-1]] +
                                 ratio.1*(levels[[i]]-levels[[i-1]])
                     }
                     ratio.n <- NA
                     y.1.n <- NA
                     z.1.n <- NA
                     if ( xy[ind.x] < isolines[[i]][[nam.x]][n.2] )
                     {
                        ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                                   (isolines[[i-1]][[nam.x]][n.1]-
                                                isolines[[i]][[nam.x]][n.2])
                        y.1.n <- isolines[[i-1]][[nam.y]][n.1] +
                                 ratio.n*(isolines[[i]][[nam.y]][n.2]-
                                                isolines[[i-1]][[nam.y]][n.1])
                        z.1.n <- levels[[i-1]] +
                                 ratio.n*(levels[[i]]-levels[[i-1]])
                     }
                     val <- utility.intpol.multiple(x  = xy[ind.y],
                                                    xs = c(y.1.1,y.1.n,y.2),
                                                    ys = c(z.1.1,z.1.n,
                                                             rep(levels[i],2)))
                     if ( ! is.na(val) )
                     {
                        z[lead.current] <- val
                        break
                     }
                  }
                  else # xy[ind.x] < min(isolines[[i-1]][[nam.x]])
                  {
                     # x coordinate of xy intersects level i but is smaller than
                     # minimum x at level i-1

                     ratio.1 <- NA
                     y.1.1 <- NA
                     z.1.1 <- NA
                     if ( xy[ind.x] > isolines[[i]][[nam.x]][1] )
                     {
                        ratio.1 <- (xy[ind.x]-isolines[[i-1]][[nam.x]][1])/
                                   (isolines[[i]][[nam.x]][1]-
                                                isolines[[i-1]][[nam.x]][1])
                        y.1.1 <- isolines[[i-1]][[nam.y]][1] +
                                 ratio.1*(isolines[[i]][[nam.y]][1]-
                                                isolines[[i-1]][[nam.y]][1])
                        z.1.1 <- levels[[i-1]] +
                                 ratio.1*(levels[[i]]-levels[[i-1]])
                     }
                     ratio.n <- NA
                     y.1.n <- NA
                     z.1.n <- NA
                     if ( xy[ind.x] > isolines[[i]][[nam.x]][n.2] )
                     {
                        ratio.n <- (isolines[[i-1]][[nam.x]][n.1]-xy[ind.x])/
                                   (isolines[[i-1]][[nam.x]][n.1]-
                                                isolines[[i]][[nam.x]][n.2])
                        y.1.n <- isolines[[i-1]][[nam.y]][n.1] +
                                 ratio.n*(isolines[[i]][[nam.y]][n.2]-
                                                isolines[[i-1]][[nam.y]][n.1])
                        z.1.n <- levels[[i-1]] +
                                 ratio.n*(levels[[i]]-levels[[i-1]])
                     }
                     val <- utility.intpol.multiple(x  = xy[ind.y],
                                                    xs = c(y.1.1,y.1.n,y.2),
                                                    ys = c(z.1.1,z.1.n,
                                                            rep(levels[i],2)))
                     if ( ! is.na(val) )
                     {
                        z[lead.current] <- val
                        break
                     }
                  }
               }
               else # not within ranges of contour lines at level i-1 and i
               {
                  x.1.1 <- isolines[[i-1]][[nam.x]][1]
                  x.2.1 <- isolines[[i]][[nam.x]][1] 
                  x.1.n <- isolines[[i-1]][[nam.x]][n.1]
                  x.2.n <- isolines[[i]][[nam.x]][n.2]
                  if ( (xy[ind.x] >= x.1.1 & xy[ind.x] <= x.2.1) |
                       (xy[ind.x] >= x.2.1 & xy[ind.x] <= x.1.1) )
                  {
                     if ( (xy[ind.x] >= x.1.n & xy[ind.x] <= x.2.n) |
                          (xy[ind.x] >= x.2.n & xy[ind.x] <= x.1.n) )
                     {
                        # x not within the ranges of isolines at lev- i-1 and i;
                        # x within the range of the bounding lines between the
                        # ends of the isolines at levels i-1 and i

                        ratio.1 <- (xy[ind.x]-x.1.1)/(x.2.1-x.1.1)
                        y.1 <- isolines[[i-1]][[nam.y]][1] + 
                               ratio.1*(isolines[[i]][[nam.y]][1]-
                                                    isolines[[i-1]][[nam.y]][1])
                        z.1 <- levels[i-1] + ratio.1*(levels[i]-levels[i-1])
                        ratio.n <- (xy[ind.x]-x.1.n)/(x.2.n-x.1.n)
                        y.n <- isolines[[i-1]][[nam.y]][n.1] + 
                               ratio.n*(isolines[[i]][[nam.y]][n.2]-
                                                  isolines[[i-1]][[nam.y]][n.1])
                        z.n <- levels[i-1] + ratio.n*(levels[i]-levels[i-1])
                        if ( (xy[ind.y] >= y.1 & xy[ind.y] <= y.n) |
                             (xy[ind.y] <= y.1 & xy[ind.y] >= y.n) )
                        {
                           z[lead.current] <- 
                                       z.1 + (xy[ind.y]-y.1)/(y.n-y.1)*(z.n-z.1)
                           break
                        }
                     }
                  }
               }
            } 
         }     
      }  
   }
   if ( is.na(z[1]) & is.na(z[2]) ) return(NA)
     
   return(mean(z,na.rm=TRUE))
}


utility.check.required <- function(u,required,num.required)
{
   res.ok <- sum(ifelse(is.na(u),0,1)) >= num.required &
             sum(ifelse(is.na(u) & required,1,0)) == 0
   return(res.ok)
}


utility.check.name <- function(name,nodes)
{
   nodes.local <- nodes
   if ( !is.list(nodes) ) nodes.local <- as.list(nodes)
   for ( i in 1:length(nodes) )
   {
      if ( name == nodes[[i]]$name ) return(FALSE)
   }
   return(TRUE)
}


utility.structure <- function(node)
{
   if ( substring(class(node),1,7) != "utility" )
   {
      warning("Node \"",node$name,"\": argument must be a subclass of utility")
      return(NA)
   }
   str <- data.frame(upper        = NA,
                     utility      = node$utility,
                     required     = node$required,
                     num.required = node$num.required,
                     color        = node$col,
                     endnode      = FALSE,
                     attributes   = NA,
                     level        = 1 + node$shift.levels,
                     endnodes     = 0,
                     offset       = 0)
   rownames(str) <- node$name
   offset <- 0
   for ( i in 1:length(node$nodes) )
   {
      if ( node$nodes[[i]]$type == "endnode" )
      {
         str.new <- data.frame(upper        = node$name,
                               utility      = node$nodes[[i]]$utility,
                               required     = node$nodes[[i]]$required,
                               num.required = NA,
                               color        = node$nodes[[i]]$col,
                               endnode      = TRUE,
                               attributes   = paste(node$nodes[[i]]$attrib,
                                                    collapse=";"),
                               level        = 2 + node$shift.levels,
                               endnodes     = 1,
                               offset       = offset)
         rownames(str.new) <- node$nodes[[i]]$name
         str[1,"endnodes"] <- str[1,"endnodes"] + 1
         offset <- offset + 1
      }
      else
      {
         str.new <- utility.structure(node$nodes[[i]])
         if ( ! is.data.frame(str.new) ) return(NA)
         str.new[1,"upper"] <- node$name
         str.new$level <- str.new$level + 1 + node$shift.levels
         str.new$offset <- str.new$offset + offset
         str[1,"endnodes"] <- str[1,"endnodes"] + str.new[1,"endnodes"]
         offset <- offset + sum(ifelse(str.new$endnode,1,0))
      }
      ind1 <- match(rownames(str.new),rownames(str))
      ind2 <- ind1[!is.na(ind1)]
      if ( length(ind2) > 0 )
      {
         cat("*** Warning: node name(s) not unique:","\n",
             paste(rownames(str)[ind2],"\n"))
         return(NA)
      }
      str <- rbind(str,str.new)
   }
   return(str)
}


utility.endnode.plot1d <- 
                   function(node,
                            col       = utility.calc.colors(),
                            gridlines = c(0.2,0.4,0.6,0.8),
                            main      = "",
                            cex.main  = 1,
                            ...)
{
   length <- 101
   x <- seq(node$range[1],node$range[2],length=length)
   u <- evaluate(node,attrib=x)
   title <- main; if ( nchar(title) == 0 ) title <- node$name
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   plot(numeric(0),numeric(0),type="l",
        xlim=node$range,ylim=c(0,1),
        xlab=node$attrib,ylab=funtype,main=title,
        xaxs="i",yaxs="i",cex.main=cex.main,...)
        
   if ( length(col)>1 & !node$utility )
   {
      # plot colored axes
      
      midpoints <- 0.5*(u[-1]+u[-length(u)])
      cols <- utility.get.colors(u,col)
      for ( i in 1:(length-1) )
      {
         lines(c(x[i],x[i+1]),c(0,0)+0.001,col=cols[i],lwd=3)
      }
      du <- 1/(length-1)
      midpoints <- seq(du,1-du,length=length-1)
      cols <- utility.get.colors(midpoints,col)
      for ( i in 1:(length-1) )
      {
         lines(c(1,1)*(node$range[1]+0.001*(node$range[2]-node$range[1])),
               c((i-1)*du,i*du),col=cols[i],lwd=3)
      }
      
      # plot grid lines:
      
      if ( ! is.na(gridlines[1]) )
      {
         for ( level in gridlines )
         {
            abline(h=level,lty="dashed")
            for ( i in 1:(length-1) )
            {
               if ( !is.na(u[i]) & !is.na(u[i+1]) )
               {
                  if ( (u[i] <= level & u[i+1] > level) |
                       (u[i] > level & u[i+1] <= level) )
                  {
                     x.level <- x[i] + (level-u[i])/(u[i+1]-u[i])*(x[i+1]-x[i])
                     lines(c(x.level,x.level),c(0,level),lty="dashed")
                  }
               }
            }
         }
      }
   }
   
   # plot value/utility function:
   
   color <- "black"
   if ( length(col) == 1 ) color <- col
   lines(x,u,lwd=2,col=color)
}


utility.endnode.plot2d <- function(node,
                                   col       = utility.calc.colors(),
                                   gridlines = c(0.2,0.4,0.6,0.8),
                                   main      = "",
                                   cex.main  = 1,
                                   ...)
{
   num.grid <- 100
   x <- node$ranges[[1]][1] + 
        ((1:num.grid)-0.5)/num.grid*(node$ranges[[1]][2]-node$ranges[[1]][1])
   y <- node$ranges[[2]][1] + 
        ((1:num.grid)-0.5)/num.grid*(node$ranges[[2]][2]-node$ranges[[2]][1])
   
   array.x <- sort(rep(x,num.grid))
   array.y <- rep(y,num.grid)
   array.xy <- cbind(array.x,array.y)
   colnames(array.xy) <- node$attrib
   
   u <- evaluate(node,as.data.frame(array.xy))
   u <- t(matrix(u,ncol=num.grid,byrow=FALSE))
   
   title <- main; if ( nchar(title) == 0 ) title <- node$name
   image(x=x,y=y,z=u,xlim=node$ranges[[1]],ylim=node$ranges[[2]],zlim=c(0,1),
         col=col,xlab=node$attrib[1],ylab=node$attrib[2],main=title,
         cex.main=cex.main)
}


utility.conversion.plot <- function(node,
                                    col       = "black",
                                    gridlines = NA,
                                    cex.main  = 1,
                                    ...)
{
   length <- 101
   x <- ((1:length)-1)/(length-1)
   u <- evaluate.cond(node,x)
   plot(numeric(0),numeric(0),type="l",
        xlim=c(0,1),ylim=c(0,1),
        xlab=paste("value(",node$nodes[[1]]$name,")",sep=""),ylab="utility",
        main=node$name,xaxs="i",yaxs="i",cex.main=cex.main)
   color <- "black"; if ( length(col) == 1 ) color <- col
   lines(x,u,lwd=2,col=color)
   lines(c(0,1),c(0,1))
   if ( length(node$x) > 0 & length(node$u) > 0 )
   {
      if ( length(node$x) == length(node$u) )
      {
         points(node$x,node$u,cex=1.5,xpd=TRUE)
      }
   }
}


utility.aggregation.plot <- function(node           = node,
                                     col            = col,
                                     gridlines      = gridlines,
                                     cex.main       = 1,
                                     ...)
{
   nodes.names <- rep(NA,length(node$nodes))
   for ( i in 1:length(node$nodes) ) nodes.names[i] <- node$nodes[[i]]$name
   if ( node$name.fun == "utility.aggregate.add" )
   {
      w <- node$par/sum(node$par)
      w.max <- max(w)
      if ( length(w) != length(nodes.names) )
      {
         warning("Node \"",node$name,"\": ",
                 "length of sub-nodes and weights not equal: ",
                 length(nodes.names)," ",length(w),sep="")
      }
      else
      {
         barplot(w,names.arg=nodes.names,ylim=c(0,1.2*w.max),
                 ylab="weight",main=node$name,cex.main=cex.main)
         text(0.5*1.3*length(w),1.1*w.max,"additive aggregation with weights:")
      }
   }
   else
   {
      if ( node$name.fun == "utility.aggregate.cobbdouglas" )
      {
         w <- node$par/sum(node$par)
         w.max <- max(w)
         if ( length(w) != length(nodes.names) )
         {
            warning("Node \"",node$name,"\": ",
                    "length of sub-nodes and weights not equal ",
                    length(nodes.names)," ",length(w),sep="")
         }
         else
         {
            barplot(w,names.arg=nodes.names,ylim=c(0,1.2*w.max),
                    ylab="weight",main=node$name,cex.main=cex.main)
            text(0.5*1.3*length(w),1.1*w.max,
                 "Cobb-Douglas aggregation with weights:")
         }
      }
      else
      {
         if ( node$name.fun == "utility.aggregate.mult" )
         {
            w <- node$par
            w.max <- max(w)
            if ( length(w) != length(nodes.names) )
            {
               warning("Node \"",node$name,"\": ",
                       "length of sub-nodes and weights not equal: ",
                       length(nodes.names)," ",length(w),sep="")
            }
            else
            {
               barplot(w,names.arg=nodes.names,ylim=c(0,1.2*w.max),
                       ylab="weight",main=node$name,cex.main=cex.main)
               text(0.5*1.3*length(w),1.1*w.max,
                    "multiplicative aggregation with weights:")
            }
         }
         else
         {
            if ( node$name.fun == "utility.aggregate.min" )
            {
               plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
                    xaxt="n",yaxt="n",main=node$name,xlab="",ylab="",
                    cex.main=cex.main)
               text(0.5,0.9,"Minimum (worst-case) aggregation of nodes:")
               for ( i in 1:length(nodes.names) )
               {
                  text(0.5,0.7*i/length(nodes.names),nodes.names[i])
               }
            }
            else
            {
               plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
                    xaxt="n",yaxt="n",main=node$name,xlab="",ylab="",
                    cex.main=cex.main)
               text(0.5,0.9,paste("aggregation with function \"",
                                  node$name.fun,"\" of nodes:",sep=""))
               for ( i in 1:length(nodes.names) )
               {
                  text(0.5,0.7*i/length(nodes.names),nodes.names[i])
               }
            }
         }
      }
   }
}


utility.plothierarchy <- 
   function(node,
            u           = NA,
            uref        = NA,
            col         = utility.calc.colors(),
            main        = main,
            cex.main    = 1,
            cex.nodes   = 1,
            cex.attrib  = 1,
            with.attrib = TRUE,
            ...)
{
   # global parameters:

   delta.x        <- 0.1
   delta.y        <- 0.1
   min.median.dev <- 0.03
   num.stripes    <- 500
   dh.rel.utility <- 0.1

   # get hierarchy structure and define positions of boxes:
         
   str <- utility.structure(node)
   if ( ! is.data.frame(str) ) return(NA)
   str$level <- str$level-min(str$level)+1  # remove indent of top node
   w <- 1/max(str$level)
   if ( with.attrib ) w <- 1/(max(str$level)+1)
   h <- 1/str$endnodes[1]
   str$x <- (str$level-0.5)*w
   str$y <- 1-(str$offset+0.5*str$endnodes)*h
   x.attrib <- max(str$level)*w + delta.y*w 

   # convert u and uref to data frames:

   u.local <- u
   if ( is.vector(u.local) ) u.local <- t(u.local)         
   u.local <- as.data.frame(u.local)
   uref.local <- uref
   if ( is.vector(uref.local) ) uref.local <- t(uref.local)         
   uref.local <- as.data.frame(uref.local)
   
   # plot indvidual plots per row if the same number of titles is provided;
   # plot quantile summary if not the same number of titles is provided and 
   # if the number of rows is > 1
   
   quant.summary <- length(main) != nrow(u.local) & nrow(u.local) > 1
   
   # find out if u and uref are available (otherwise plot required/not required shading)

   u.available <- FALSE
   if ( nrow(u.local)>1 | ncol(u.local)>1 | !is.na(u.local[1,1]) )
   {
      u.available <- TRUE
   }
   uref.available <- FALSE 
   ind.uref.local <- rep(1,nrow(u.local))
   if ( nrow(uref.local)>1 | ncol(uref.local)>1 | !is.na(uref.local[1,1]) )
   {
      uref.available <- TRUE
      if ( !quant.summary ) # number of rows must be unity or equal to nrow(u)
      {
         if ( nrow(uref.local) == nrow(u.local) )
         {
            ind.uref.local <- 1:nrow(u.local)
         }
         else
         {
            if ( nrow(uref.local) != 1 ) uref.available <- FALSE
         }
      }
   }
   
   # loop over rows of utilities/values:

   num.plots <- nrow(u.local)
   if ( !u.available | quant.summary ) num.plots <- 1
   for ( k in 1:num.plots )
   {
      # set-up plot frame:
         
      par.def <- par(no.readonly=TRUE)
      par(mar=c(0,0,0,0))
      plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
           xaxt="n",yaxt="n",xlab="",ylab="",cex.main=cex.main)
           
      # write title
      
      title <- main[1]
      if ( length(main) == nrow(u.local) ) title <- main[k]
      text(0,1-0.5*h,title,adj=c(0,0.5),cex=cex.main,...)
      
      # draw color code legend:
      
      if ( u.available )
      {
         x.l <- delta.x*w
         x.r <- (1-delta.x)*w
         y   <- 0.8*h
         num.col <- 100
         v <- (1:num.col - 0.5)/num.col
         colors <- utility.get.colors(v,col)
         for ( i in 1:num.col ) 
         {
            lines(x.l+(x.r-x.l)/num.col*c(i-1,i),c(y,y),col=colors[i],lwd=3)
         }
         text(x.l,y,"0",pos=1,cex=cex.nodes)
         text(x.r,y,"1",pos=1,cex=cex.nodes)
      }
      
      for ( i in 1:nrow(str) )
      {
         # calculate box edge coordinates:
            
         x.box.l <- str$x[i] - (0.5-delta.x)*w
         x.box.r <- str$x[i] + (0.5-delta.x)*w
         y.box.b <- str$y[i] - (0.5-delta.y)*h
         y.box.t <- str$y[i] + (0.5-delta.y)*h

         # plot background color:
            
         if ( !u.available ) # plot required/not required nodes in differnt grey
         {
            if ( str$required[i] ) color <- grey(0.7)
            else                   color <- grey(0.9)
            polygon(x   = c(x.box.l,x.box.r,x.box.r,x.box.l,x.box.l),
                    y   = c(y.box.b,y.box.b,y.box.t,y.box.t,y.box.b),
                    col = color)
         }
         else
         {
            if ( !quant.summary ) # plot hierarchy for each row of u
            {
               if ( str$utility[i] )  # plot horizontal lines for utility node
               {
                  val <- u.local[k,rownames(str)[i]]
                  if ( ! is.na(val) )
                  {
                    y.t <- y.box.t
                    if ( uref.available ) y.t <- 0.5*(y.box.b+y.box.t)
                    lines((x.box.l+val*(x.box.r-x.box.l)*c(1,1)),
                           c(y.box.b,y.t),lwd=1.5)
                  }
                  if ( uref.available )
                  {
                     val <- uref.local[ind.uref.local[k],rownames(str)[i]]
                     if ( ! is.na(val) )
                     {
                        y.b <- 0.5*(y.box.b+y.box.t)
                        lines((x.box.l+val*(x.box.r-x.box.l)*c(1,1)),
                               c(y.b,y.box.t),lwd=1.5)
                     }
                  }
               }
               else   # plot colored boxes for value nodes
               {
                  color <- "white"
                  val <- u.local[k,rownames(str)[i]]
                  if ( ! is.na(val) )
                  {
                     color <- utility.get.colors(val,col)
                  }
                  y.t <- y.box.t
                  if ( uref.available ) y.t <- 0.5*(y.box.b+y.box.t)
                  polygon(x   = c(x.box.l,x.box.r,x.box.r,x.box.l,x.box.l),
                          y   = c(y.box.b,y.box.b,y.t    ,y.t    ,y.box.b),
                          col = color,border=NA)
                  # plot black value line:
                  lines((x.box.l+val*(x.box.r-x.box.l)*c(1,1)),
                        c(y.box.b,y.t),lwd=1.5)
                  if ( uref.available )
                  {
                     color <- "white"
                     val <- uref.local[ind.uref.local[k],rownames(str)[i]]
                     if ( ! is.na(val) )
                     {
                        color <- utility.get.colors(val,col)
                     }
                     y.b <- 0.5*(y.box.b+y.box.t)
                     polygon(x   = c(x.box.l,x.box.r,x.box.r,x.box.l,x.box.l),
                             y   = c(y.b    ,y.b    ,y.box.t,y.box.t,y.b),
                             col = color,border=NA)
                     # plot black value line:
                     lines((x.box.l+val*(x.box.r-x.box.l)*c(1,1)),
                           c(y.b,y.box.t),lwd=1.5)
                  }
               }
            }
            else # plot quantile summary of v or expected u
            {
               if ( str$utility[i] )  # plot vertical line for expected utility
               {
                  u.exp <- NA
                  column <- match(rownames(str)[i],colnames(u.local))
                  if ( !is.na(column) )
                  {
                     u.exp <- mean(u.local[,column],na.rm=TRUE)
                  }
                  uref.exp <- NA
                  if ( uref.available )
                  {
                     column <- match(rownames(str)[i],colnames(uref.local))
                     if ( !is.na(column) )
                     {
                        uref.exp <- mean(uref.local[,column],na.rm=TRUE)
                     }
                  }
                  if ( ! ( is.na(u.exp) | is.na(uref.exp) ) ) # illustrate better
                  {                                           # alternative
                     color     <- "lightgreen"
                     color.ref <- "tomato"
                     if ( uref.exp > u.exp )
                     {
                        color     <- "tomato"
                        color.ref <- "ligthgreen"
                     }
                     y.t <- y.box.t
                     if ( uref.available ) y.t <- 0.5*(y.box.b+y.box.t)
                     polygon(x   = c(x.box.l,x.box.r,x.box.r,x.box.l,x.box.l),
                             y   = c(y.box.b,y.box.b,y.t    ,y.t    ,y.box.b),
                                   col = color,border=NA)
                     y.b <- 0.5*(y.box.b+y.box.t)
                     polygon(x   = c(x.box.l,x.box.r,x.box.r,x.box.l,x.box.l),
                             y   = c(y.b    ,y.b    ,y.box.t,y.box.t,y.b),
                             col = color.ref,border=NA)
                  }
                  if ( ! is.na(u.exp) )
                  {
                     y.t <- y.box.t
                     if ( uref.available ) y.t <- 0.5*(y.box.b+y.box.t)
                     lines((x.box.l+u.exp*(x.box.r-x.box.l)*c(1,1)),
                           c(y.box.b,y.t),lwd=1.5)
                  }
                  if ( uref.available )
                  {
                     if ( ! is.na(uref.exp) )
                     {
                       y.b <- 0.5*(y.box.b+y.box.t)
                       lines((x.box.l+uref.exp*(x.box.r-x.box.l)*c(1,1)),
                             c(y.b,y.box.t),lwd=1.5)
                     }
                  }
               }
               else  # plot quantile summary for value nodes
               {
                  column <- match(rownames(str)[i],colnames(u.local))
                  if ( !is.na(column) )
                  {
                     u.vals <- u.local[,column]
                     if ( sum(!is.na(u.vals)) > 0 )
                     {
                        u.vals <- u.vals[!is.na(u.vals)]
                        u.quant <- quantile(u.vals,probs=c(0.05,0.5,0.95))
                        if ( u.quant[2]-u.quant[1] < min.median.dev )
                        {
                           u.quant[1] <- max(0,u.quant[2]-min.median.dev )
                        }
                        if ( u.quant[3]-u.quant[2] < min.median.dev )
                        {
                           u.quant[3] <- min(1,u.quant[2]+min.median.dev )
                        }
                        y.t <- y.box.t
                        if ( uref.available ) y.t <- 0.5*(y.box.b+y.box.t)

                        # plot colored 90% credibility interval:
                        for ( j in floor(num.stripes*u.quant[1]):ceiling(num.stripes*u.quant[3]) )
                        {
                           lines((x.box.l+j/num.stripes*(x.box.r-x.box.l)*c(1,1)),
                                 c(y.box.b,y.t),
                                 col=utility.get.colors(j/num.stripes,col))
                        }
                     
                        # plot black median line:
                        lines((x.box.l+u.quant[2]*(x.box.r-x.box.l)*c(1,1)),
                              c(y.box.b,y.t),lwd=1.5)
                     }
                  }
                  if ( uref.available )
                  {
                     column <- match(rownames(str)[i],colnames(uref.local))
                     if ( !is.na(column) )
                     {
                        uref.vals <- uref.local[,column]
                        if ( sum(!is.na(uref.vals)) > 0 )
                        {
                           uref.vals <- uref.vals[!is.na(uref.vals)]
                           uref.quant <- quantile(uref.vals,probs=c(0.05,0.5,0.95))
                           if ( uref.quant[2]-uref.quant[1] < min.median.dev )
                           {
                              uref.quant[1] <- max(0,uref.quant[2]-min.median.dev )
                           }
                           if ( uref.quant[3]-uref.quant[2] < min.median.dev )
                           {
                              uref.quant[3] <- min(1,uref.quant[2]+min.median.dev )
                           }
                           y.b <- 0.5*(y.box.b+y.box.t)

                           # plot colored 90% credibility interval:
                           for ( j in floor(num.stripes*uref.quant[1]):ceiling(num.stripes*uref.quant[3]) )
                           {
                              lines((x.box.l+j/num.stripes*(x.box.r-x.box.l)*c(1,1)),
                                    c(y.b,y.box.t),
                                    col=utility.get.colors(j/num.stripes,col))
                           }
                     
                           # plot black median line:
                           lines((x.box.l+uref.quant[2]*(x.box.r-x.box.l)*c(1,1)),
                                 c(y.b,y.box.t),lwd=1.5)
                        }
                     }
                  }
               }
            }
         }
                        
         # plot bounding box:

         lines(x   = c(x.box.l,x.box.r,x.box.r,x.box.l,x.box.l),
               y   = c(y.box.b,y.box.b,y.box.t,y.box.t,y.box.b),
               col = as.character(str$color[i]))
         if ( str$utility[i] )
         {
            dh <- dh.rel.utility*(y.box.t-y.box.b)
            lines(c(x.box.l,x.box.r),(y.box.b+dh)*c(1,1))
            lines(c(x.box.l,x.box.r),(y.box.t-dh)*c(1,1))
         }
                  
         # write text into box:
            
         text(str$x[i],str$y[i],rownames(str)[i],cex=cex.nodes,...)

         # plot connecting lines:
                           
         upper <- str$upper[i]
         if ( ! is.na(upper) )
         {
            x.line.l <- str[upper,"x"] + (0.5-delta.x)*w
            x.line.r <- str$x[i] - (0.5-delta.x)*w
            x.line.v <- str[upper,"x"] + 0.5*w
            y.line.l <- str[upper,"y"]
            y.line.r <- str$y[i]
            lines(x = c(x.line.l,x.line.v,x.line.v,x.line.r), 
                  y = c(y.line.l,y.line.l,y.line.r,y.line.r))
         }
            
         # write attribute names:
                 
         if ( with.attrib )
         {
            if ( str$endnode[i] )
            {
               attributes <- strsplit(str$attributes[i],split=";")[[1]]
               n <- length(attributes)
               for ( j in 1:n )
               {
                  y.attrib <- str$y[i] +  (0.5 - (j-0.5)/n)*(1-delta.y)*h
                  text(x.attrib,y.attrib,attributes[j],pos=4,cex=cex.attrib,...)
                  lines(c(x.box.r,x.attrib),c(y.attrib,y.attrib),lty="dotted")
               }
            }
         }
      } # end for i
      par(par.def)
   } # end for k
}


utility.plottable <- 
   function(u          = NA,
            nodes      = NA,
            col        = utility.calc.colors(),
            main       = "",
            labels     = NA,
            cex.main   = 1,
            cex.nodes  = 1,
            f.reaches  = 0.2,
            f.nodes    = 0.2,
            ...)
{
   # global parameters:

   delta.x    <- 0.2
   delta.y    <- 0.2
   delta.main <- 0.05

   # initializations:
     
   if ( length(dim(u)) != 2 ) return
   ind.reaches <- 1:nrow(u)
   ind.nodes <- 1:ncol(u)
   if ( !is.na(nodes[1]) )
   {
      ind.nodes <- match(nodes,colnames(u))
   }
   
   # set-up plotting parameters and plot frame:

   dx <- (1-f.reaches)/length(ind.nodes)
   dy <- (1-f.nodes)/length(ind.reaches)
   x <- f.reaches+(1:length(ind.nodes)-0.5)*dx
   y <- 1-f.nodes-(1:length(ind.reaches)-0.5)*dy
   if ( nchar(main[1]) > 0 )
   {
      y  <- (1-delta.main)*y
      dy <- (1-delta.main)*dy
   }
   par.def <- par(no.readonly=TRUE)
   par(mar=c(0,0,0,0))
   plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
        xaxt="n",yaxt="n",xlab="",ylab="")
   
   # write names of reaches:

   if ( nchar(main[1]) > 0 ) text(x=0.5,y=1-0.5*delta.main,label=main[1],cex=cex.main)
   for ( i in 1:length(ind.reaches) ) 
   {
      if ( !is.na(ind.reaches[i]) )
      {
         if ( is.na(labels[1]) )
         {
            text(x=0,y=y[i],label=rownames(u)[ind.reaches[i]],adj=c(0,0.5),cex=cex.nodes)
         }
         else
         {
            text(x=0,y=y[i],label=labels[ind.reaches[i]],adj=c(0,0.5),cex=cex.nodes)
         }
      }
   }
   
   # write and color values:
   
   for ( i in 1:length(ind.reaches) )
   {
      for ( j in 1:length(ind.nodes) )
      {
         if ( !is.na(ind.reaches[i]) & !is.na(ind.nodes[j]) )
         {
            color <- "white"
            val.str <- ""
            val <- u[ind.reaches[i],ind.nodes[j]]
            if ( ! is.na(val) )
            {
               color <- utility.get.colors(val,col)
               val.str <- paste(round(val,2))
               if ( nchar(val.str) > 1 & substring(val.str,1,1) == "0" )
               {
                  val.str <- substring(val.str,2)
                  if ( nchar(val.str) == 2 ) val.str <- paste(val.str,"0",sep="")
               }
            }
            polygon(x   = x[j]+0.5*(1-delta.x)*dx*c(-1,1,1,-1,-1),
                    y   = y[i]+0.5*(1-delta.y)*dy*c(1,1,-1,-1,1),
                    col = color)
            text(x=x[j],y=y[i],val.str,cex=cex.nodes)
         }
      }
   }

   
   # write names of nodes:
   
   par(srt=90)
   for ( j in 1:length(ind.nodes) ) 
   {
      if ( !is.na(ind.nodes[j]) )
      {
         text(x=x[j],y=1-f.nodes,
              label=colnames(u)[ind.nodes[j]],
              adj=c(0,0.5),cex=cex.nodes)
      }
   }
   
   # reset plotting parameters:
   
   par(par.def)
}


utility.plot <- function(node,
                         u           = NA,
                         uref        = NA,
                         type        = c("hierarchy","table","node","nodes"),
                         nodes       = NA,
                         col         = utility.calc.colors(),
                         gridlines   = c(0.2,0.4,0.6,0.8),
                         main        = "",
                         labels      = NA,
                         cex.main    = 1,
                         cex.nodes   = 1,
                         cex.attrib  = 1,
                         f.reaches   = 0.2,
                         f.nodes     = 0.2,
                         with.attrib = TRUE,
                         ...)
{
   if ( type[1] == "nodes" | type[1] == "node" )
   {
      # plot current node:
      
      if ( is.na(nodes[1]) | ! is.na(match(node$name,nodes)) )
      {
         if ( substring(class(node),1,18) == "utility.conversion" )
         {
            utility.conversion.plot(node       = node,
                                    col        = col,
                                    gridlines  = gridlines,
                                    cex.main   = cex.main,
                                    cex.nodes  = cex.nodes,
                                    cex.attrib = cex.attrib,
                                    ...)
         }
         else
         {
            if ( substring(class(node),1,19) == "utility.aggregation" )
            {
              utility.aggregation.plot(node       = node,
                                        col        = col,
                                        gridlines  = gridlines,
                                        cex.main   = cex.main,
                                        cex.nodes  = cex.nodes,
                                        cex.attrib = cex.attrib,
                                        ...)
            }
            else
            {
               if ( node$type == "endnode" )
               {
                  if ( class(node) == "utility.endnode.cond" )
                  {
                    plot(node$nodes[[i]],
                          par       = NA,
                          col       = col,
                          gridlines = gridlines,
                          cex.main  = cex.main,
                          nodes     = nodes,
                          ...)
                  }
                  else
                  {
                     plot(node$nodes[[i]],
                          par       = NA,
                          col       = col,
                          gridlines = gridlines,
                          cex.main  = cex.main,
                          ...)
                  }
               }
               else
               {
                  # unknown node type; not plotted
               }
            }
         }
      }
      
      # plot other nodes:
      
      if ( type == "nodes" )
      {
         if ( length(node$nodes) > 0 )
         {
         for ( i in 1:length(node$nodes) )
         {
            # initiate plot of subnodes:

            if ( node$nodes[[i]]$type == "endnode" )
            {
              if ( class(node$nodes[[i]]) == "utility.endnode.cond" )
              {
                plot(node$nodes[[i]],
                     par       = NA,
                     col       = col,
                     gridlines = gridlines,
                     cex.main  = cex.main,
                     nodes     = nodes,
                     ...)
              }
              else
              {
                 if ( is.na(nodes[1]) | ! is.na(match(node$nodes[[i]]$name,nodes)) )
                 {
                    plot(node$nodes[[i]],
                         par       = NA,
                         col       = col,
                         gridlines = gridlines,
                         cex.main  = cex.main,
                         ...)
                  }
               }
            }
            else
            {
               plot(node$nodes[[i]],
                    u          = u,
                    par        = NA,
                    type       = type,
                    nodes      = nodes,
                    col        = col,
                    gridlines  = gridlines,
                    cex.main   = cex.main,
                    cex.nodes  = cex.nodes,
                    cex.attrib = cex.attrib,
                    ...)
            }
         }
         }
      }
   }
   else
   {
      if ( type[1] == "hierarchy" )
      {
         if ( is.na(nodes[1]) | ! is.na(match(node$name,nodes)) )
         {
            utility.plothierarchy(node        = node,
                                  u           = u,
                                  uref        = uref,
                                  col         = col,
                                  main        = main,
                                  cex.main    = cex.main,
                                  cex.nodes   = cex.nodes,
                                  cex.attrib  = cex.attrib,
                                  with.attrib = with.attrib,
                                  ...)
         }
         if ( ! is.na(nodes[1]) )
         {
            if ( node$type != "endnode" )
            {
               for ( i in 1:length(node$nodes) )
               {
                  utility.plot(node$nodes[[i]],
                               u          = u,
                               uref       = uref,
                               type       = type,
                               nodes      = nodes,
                               col        = col,
                               gridlines  = gridlines,
                               main       = main,
                               cex.main   = cex.main,
                               cex.nodes  = cex.nodes,
                               cex.attrib = cex.attrib,
                               ...)
               } 
            }
         }
      }
      else
      {
         if ( type[1] == "table" )
         {
            if ( length(dim(u)) == 2 )
            {
               utility.plottable(u          = u,
                                 nodes      = nodes,
                                 col        = col,
                                 main       = main,
                                 labels     = labels,
                                 cex.main   = cex.main,
                                 cex.nodes  = cex.nodes,
                                 f.reaches  = f.reaches,
                                 f.nodes    = f.nodes,
                                 ...)
            }
         }
         else
         {
            cat("unknown plot type:",type[1],"\n")
         }
      }
   }
}


# ==============================================================================
# endnode for single-attribute interpolation: 
# class "utility.endnode.intpol1d"
# ==============================================================================


# constructor:
# ------------

utility.endnode.intpol1d.create <- function(name.node,    # character(1)
                                            name.attrib,  # character(1)
                                            range,        # numeric(2)
                                            x,            # numeric(n)
                                            u,            # numeric(n)
                                            names.x     = rep(NA,length(x)),
                                            names.u     = rep(NA,length(u)),
                                            utility     = TRUE,
                                            required    = FALSE,
                                            col         = "black",
                                            shift.levels = 0)
{
   # consistency checks:

   check.ok <- T   
   if ( length(x) != length(u) )
   {
      cat("*** Warning: x and u of different length:",
              length(x),length(u))
      check.ok <- F
   }
   if ( length(names.x) != length(names.u) )
   {
      cat("*** Warning: names.x and names.u of different length:",
          length(names.x),length(names.u),"\n")
      check.ok <- F
   }
   if ( length(x) != length(names.x) )
   {
      cat("*** Warning: x and names.x of different length:",
          length(x),length(names.x),"\n")
      check.ok <- F
   }
   if ( range[1] >= range[2] )
   {
      cat("*** Warning: Minimum of range not smaller than maximum:",
          range[1],range[2],"\n")
      check.ok <- F
   }
   if ( sum(x[-1]-x[-length(x)] > 0) != length(x)-1 &
        sum(x[-1]-x[-length(x)] < 0) != length(x)-1 )
   {
      cat("*** Warning: x values in interpolation node must either be","\n",
          "strictly increasing or strictly decreasing","\n")
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
   node$name        <- name.node
   node$description <- "utility/value 1d interpolation end node" 
   node$type        <- "endnode"
   node$attrib      <- name.attrib
   node$range       <- range
   node$x           <- x
   node$u           <- u
   node$names.x     <- names.x
   node$names.u     <- names.u
   node$required    <- required
   node$utility     <- utility
   node$col         <- col
   node$shift.levels <- shift.levels
   class(node)      <- "utility.endnode.intpol1d" 
   
   # print and return class

   cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
   return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.intpol1d <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update adequate values in interpolation list:
      
   n <- node
   for ( i in 1:length(n$x) )
   {
      if ( ! is.na(n$names.x[i]) )
      {
         ind <- which(n$names.x[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter \"",
                    names(par)[ind[1]],"\"",sep="")
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$x[i] <- par[ind]
         }
      } 
      if ( ! is.na(n$names.u[i]) )
      {
         ind <- which(n$names.u[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter",
                    names(par)[ind[1]])
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$u[i] <- par[ind]
         }
      } 
   }
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.intpol1d <- function(x,
                                              attrib,   # data.frame, numeric
                                              par = NA,
                                              ...)
{
   node <- x
   
   # update parameters:

   n <- updatepar(node,par)
   
   # extract attributes:
   
   if ( is.data.frame(attrib) | is.matrix(attrib) )
   {
      if ( length(which(colnames(attrib)==n$attrib)) != 1 )
      {
         warning("Node \"",node$name,"\": attribute \"",n$attrib,"\" not found",sep="")
         return(rep(NA,nrow(attrib)))
      }
      a <- attrib[,n$attrib]
   }
   else
   {
      if ( ! is.vector(attrib) )
      {
         warning("Node \"",node$name,"\": unknown format of attribute \"",n$attrib,"\"",sep="")
         return(NA)
      }
      if ( length(names(attrib)) == 0 )
      {
         a <- attrib
      }
      else
      {
         ind <- which(names(attrib)==n$attrib)
         if ( length(ind) != 1 )
         {
            if ( length(ind) > 1)
            {
               warning("Node \"",node$name,"\": multiple occurrences of attribute \"",
                       n$attrib,"\"",sep="")
            }
            else
            {
               warning("Node \"",node$name,"\": attribute \"",n$attrib,"\" not found",sep="")
            }
            return(NA)
         }
         a <- attrib[ind]
      }
   }
   
   # evaluate results:
   
   u <- approx(x=n$x,y=n$u,xout=a)$y
   ind.out.of.range <- (a < n$range[1]) | (a > n$range[2])
   u <- ifelse(ind.out.of.range,NA,u)
   if ( sum(ind.out.of.range,na.rm=T) > 0 )
   {
      ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
      warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib,"\" out of range: ",
              paste(a[ind.not.na],collapse=","),sep="")
   }
      
   # return results:
   
   return(u)
}


# print:
# ------

print.utility.endnode.intpol1d <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}

# summary:
# --------

summary.utility.endnode.intpol1d <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   cat("attribute:      ",node$attrib,"\n")
   cat("attribute range:",node$range[1],"-",node$range[2],"\n")
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   cat("function type:  ",funtype,"\n")
   cat("required:       ",node$required,"\n")
   cat("data pairs:","\n")
   names.x <- ifelse(is.na(node$names.x),"",node$names.x)
   names.u <- ifelse(is.na(node$names.u),"",node$names.u)
   print(data.frame(names.x=names.x,x=node$x,u=node$u,names.u=names.u))
}


# plot:
# -----

plot.utility.endnode.intpol1d <- 
                   function(x,
                            par       = NA,
                            col       = utility.calc.colors(),
                            gridlines = c(0.2,0.4,0.6,0.8),
                            main      = "",
                            cex.main  = 1,
                            ...)
{
   node <- x
   n <- updatepar(node,par)
   utility.endnode.plot1d(node      = n,
                          col       = col,
                          gridlines = gridlines,
                          main      = main,
                          cex.main  = cex.main,
                          ...)
   points(n$x,n$u,cex=1.5,xpd=TRUE) 
}


# ==============================================================================
# endnode for interpolation based on isolines of two attributes: 
# class "utility.endnode.intpol2d"
# ==============================================================================


# constructor:
# ------------

utility.endnode.intpol2d.create <- function(name.node,   # character(1)
                                            name.attrib, # character(2)
                                            ranges,      # list(2) of numeric(2)
                                            isolines,    # list(n) of list of
                                                         # x, y, and, optionally
                                                         # names.x, names.y
                                            u,           # numeric(n)
                                            names.u     = rep(NA,length(u)),
                                            lead        = 0,
                                            utility     = TRUE,
                                            required    = FALSE,
                                            col         = "black",
                                            shift.levels = 0)
{
   # consistency checks:

   check.ok <- T 
   if ( length(name.attrib) != 2 )
   {
      cat("*** Warning: name.attrib must be of length 2","\n")
      check.ok <- F
   }
   if ( length(ranges) != 2 )
   {
      cat("*** Warning: ranges must be a list of two ranges","\n")
      check.ok <- F
   }
   else
   {
      if ( length(ranges[[1]]) != 2 )
      {
         cat("*** Warning: ranges[[1]] must contain two elements","\n")
         check.ok <- F
      }
      else
      {
         if ( ranges[[1]][1] >= ranges[[1]][2] )
         {
            cat("*** Warning: Minimum of range not smaller than maximum:",
                ranges[[1]][1],ranges[[1]][2],"\n")
            check.ok <- F
         }
      }
      if ( length(ranges[[2]]) != 2 )
      {
         cat("*** Warning: ranges[[2]] must contain two elements","\n")
         check.ok <- F
      }
      else
      {
         if ( ranges[[2]][1] >= ranges[[2]][2] )
         {
            cat("*** Warning: Minimum of range not smaller than maximum:",
                ranges[[2]][1],ranges[[2]][2],"\n")
            check.ok <- F
         }
      }
   } 
   if ( length(isolines) < 2 )
   {
      cat("*** Warning: at least two isolines are required","\n")
      check.ok <- F
   } 
   if ( length(isolines) != length(u) )
   {
      cat("*** Warning: isolines and u are of different length:",
          length(isolines),length(u),"\n")
      check.ok <- F
   }
   for ( i in 1:length(isolines) )
   {
      len.x <- length(isolines[[i]]$x) 
      if ( len.x < 2 )
      {
         cat("*** Warning: element x of isoline[[",i,"]] ",
             "must be of length > 1","\n",sep="")
         check.ok <- F
      }
      if ( len.x != length(isolines[[i]]$y)  )
      {
         cat("*** Warning: x and y in isoline[[",i,"]] ",
             "have different lengths:",
             len.x," ",length(isolines[[i]]$y),"\n",
             sep="")
         check.ok <- F
      }
      if ( length(isolines[[i]]$names.x) == 0 ) isolines[[i]]$names.x <- rep(NA,len.x) 
      if ( len.x != length(isolines[[i]]$names.x) ) 
      {
         cat("*** Warning: x and names.x in isoline[[",i,"]] ",
             "have different lengths:",
             len.x," ",length(isolines[[i]]$names.x),"\n",
             sep="")
         check.ok <- F
      }
      if ( length(isolines[[i]]$names.y) == 0 ) isolines[[i]]$names.y <- rep(NA,len.x) 
      if ( len.x != length(isolines[[i]]$names.y) ) 
      {
         cat("*** Warning: y and names.y in isoline[[",i,"]] ",
             "have different lengths:",
             len.x," ",length(isolines[[i]]$names.y),"\n",
             sep="")
         check.ok <- F
      }      
   }
   if ( ! check.ok )
   {
      cat("*** Warning: node \"",name.node,"\" could not be constructed","\n",
          sep="")
      return(NA)
   }

   # construct class:
   
   node <- list()
   node$name        <- name.node
   node$description <- "utility/value 2d interpolation end node"
   node$type        <- "endnode"
   node$attrib      <- name.attrib
   node$ranges      <- ranges
   node$isolines    <- isolines
   node$u           <- u
   node$names.u     <- names.u
   node$lead        <- lead
   node$required    <- required
   node$utility     <- utility
   node$col         <- col
   node$shift.levels <- shift.levels
   class(node)      <- "utility.endnode.intpol2d" 
   
   # print and return class

   cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
   return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.intpol2d <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update adequate values in interpolation list:
      
   n <- node
   for ( i in 1:length(n$u) )
   {
      if ( ! is.na(n$names.u[i]) )
      {
         ind <- which(n$names.u[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter",
                    names(par)[ind[1]],sep="")
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$u[i] <- par[ind]
         }
      } 
      for ( j in 1:length(n$isolines[[i]]$x) )
      {
         if ( ! is.na(n$isolines[[i]]$names.x[j]) )
         {
            ind <- which(n$isolines[[i]]$names.x[j] == names(par) )
            if ( length(ind) > 1 )
            {
               warning("Node \"",node$name,"\": multiple occurrences of parameter",
                       names(par)[ind[1]],sep="")
               ind <- ind[1]
            }
            if ( length(ind) == 1 )
            {
               n$isolines[[i]]$x[j] <- par[ind]
            }
         }
         if ( ! is.na(n$isolines[[i]]$names.y[j]) )
         {
            ind <- which(n$isolines[[i]]$names.y[j] == names(par) )
            if ( length(ind) > 1 )
            {
              warning("Node \"",node$name,"\": multiple occurrences of parameter",
                      names(par)[ind[1]],sep="")
               ind <- ind[1]
            }
            if ( length(ind) == 1 )
            {
               n$isolines[[i]]$y[j] <- par[ind]
            }
         }
      } 
   }
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.intpol2d <- function(x,
                                              attrib,   # data.frame, numeric
                                              par = NA,
                                              ...)
{
   node <- x
   
   # update parameters:

   n <- updatepar(node,par)
   
   # extract attributes:
   
   if ( is.data.frame(attrib) | is.matrix(attrib) )
   {
      ind <- match(n$attrib,colnames(attrib))
      if ( sum(ifelse(is.na(ind),1,0)) > 0 )
      {
         warning("Node \"",node$name,"\": attribute(s) \"",
                 paste(n$attrib[is.na(ind)],collapse=","),"\" not found",sep="")
         return(rep(NA,nrow(attrib)))
      }
      a <- attrib[,ind]
   }
   else
   {
      if ( ! is.vector(attrib) )
      {
         warning("Node \"",node$name,"\": unknown format of attribute(s) \"",n$attrib,"\"",sep="")
         return(NA)
      }
      if ( length(names(attrib)) == 0 )
      {
         if ( length(attrib) == 2 )
         a <- as.matrix(attrib,nrow=1)
      }
      else
      {
         ind <- match(n$attrib,names(attrib))
         if ( sum(ifelse(is.na(ind),1,0)) > 0 )
         {
           warning("Node \"",node$name,"\": attribute(s) \"",
                   paste(n$attrib[is.na(ind)],collapse=","),"\" not found",sep="")
           return(rep(NA,nrow(attrib)))
         }
         a <- as.matrix(attrib[ind],nrow=1)
      }
   }
   
   # evaluate results:
   
   ind <- order(n$u)
   u <- utility.intpol2d(xy=a,isolines=n$isolines[ind],
                         levels=n$u[ind],lead=n$lead)

   ind.out.of.range <- (a[,1]<n$range[[1]][1])|(a[,1]>n$range[[1]][2])
   u <- ifelse(ind.out.of.range,NA,u)
   if ( sum(ind.out.of.range,na.rm=T) > 0 )
   {
     ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
     warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib[1],"\" out of range: ",
             paste(a[ind.not.na,1],collapse=","),sep="")
   }
   
   ind.out.of.range <- (a[,2]<n$range[[2]][1])|(a[,2]>n$range[[2]][2])
   u <- ifelse(ind.out.of.range,NA,u)
   if ( sum(ind.out.of.range,na.rm=T) > 0 )
   {
     ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
     warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib[2],"\" out of range: ",
             paste(a[ind.not.na,2],collapse=","),sep="")
   }
   
   # return results:
   
   return(u)
}


# print:
# -----

print.utility.endnode.intpol2d <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.intpol2d <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   cat("attributes:      ",paste(node$attrib,collapse=" , "),"\n")
   cat("attribute ranges:",node$range[[1]][1],"-",node$range[[1]][2],
       ",",node$range[[2]][1],"-",node$range[[2]][2],"\n")
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   cat("function type:   ",funtype,"\n")
   cat("required:        ",node$required,"\n")
   cat("isolines:","\n")
   for ( i in 1:length(node$u) )
   {
      name.u <- ""
      if ( !is.na(node$names.u[i]) ) 
      {
         name.u <- paste(":",node$names.u[i])
      }
      cat("u:",node$u[i],"  ",name.u,"\n")
      names.x <- rep("",length(node$isolines[[i]]$x))
      if ( length(node$isolines[[i]]$names.x) > 0 )
      {    
         names.x <- ifelse(is.na(node$isolines[[i]]$names.x),
                           "",node$isolines[[i]]$names.x)
      }
      names.y <- rep("",length(node$isolines[[i]]$y))
      if ( length(node$isolines[[i]]$names.y) > 0 )
      {    
         names.y <- ifelse(is.na(node$isolines[[i]]$names.y),
                           "",node$isolines[[i]]$names.y)
      }
      print(data.frame(names.x=names.x,
                       x=node$isolines[[i]]$x,
                       y=node$isolines[[i]]$y,
                       names.y=names.y))
   }
}


# plot:
# -----

plot.utility.endnode.intpol2d <- 
                   function(x,
                            par       = NA,
                            col       = utility.calc.colors(),
                            gridlines = c(0.2,0.4,0.6,0.8),
                            main      = "",
                            cex.main  = 1,
                            ...)
{
   node <- x
   n <- updatepar(node,par)
   utility.endnode.plot2d(node      = n,
                          col       = col,
                          gridlines = gridlines,
                          main      = main,
                          cex.main  = cex.main,
                          ...)
   ind <- order(n$u)
   levels <- n$u[ind]
   isolines <- n$isolines[ind]
   for ( i in 1:length(levels) )
   {
      lines(isolines[[i]],...)
      if ( i > 1 )
      {
         lines(c(isolines[[i-1]]$x[1],isolines[[i]]$x[1]),
               c(isolines[[i-1]]$y[1],isolines[[i]]$y[1]),
               ...)
         lines(c(isolines[[i-1]]$x[length(isolines[[i-1]]$x)],
                 isolines[[i]]$x[length(isolines[[i]]$x)]),
               c(isolines[[i-1]]$y[length(isolines[[i-1]]$y)],
                 isolines[[i]]$y[length(isolines[[i]]$x)]),
                 ...)
      }
   }
}


# ==============================================================================
# endnode for 1d (single attribute) parametric function: 
# class "utility.endnode.parfun1d"
# ==============================================================================


# constructor:
# ------------

utility.endnode.parfun1d.create <- function(name.node,    # character(1)
                                            name.attrib,  # character(1)
                                            range,        # numeric(2)
                                            name.fun,     # name of f(a,par)
                                            par,          # numeric(n)
                                            names.par   = rep(NA,length(par)),
                                            utility     = TRUE,
                                            required    = FALSE,
                                            col         = "black",
                                            shift.levels = 0)
{
   # consistency checks:

   check.ok <- T   
   if ( length(par) != length(names.par) )
   {
      cat("*** Warning: par and names.par of different length:",
          length(par),length(names.par),"\n")
      check.ok <- F
   }
   if ( range[1] >= range[2] )
   {
      cat("*** Warning: Minimum of range not smaller than maximum:",
          range[1],range[2],"\n")
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
   node$name        <- name.node
   node$description <- "utility/value 1d parametric function end node"
   node$type        <- "endnode"
   node$attrib      <- name.attrib
   node$range       <- range
   node$name.fun    <- name.fun
   node$par         <- par
   node$names.par   <- names.par
   node$required    <- required
   node$utility     <- utility
   node$col         <- col
   node$shift.levels <- shift.levels
   class(node)      <- "utility.endnode.parfun1d" 
   
   # print and return class
   
   cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
   return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.parfun1d <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update adequate values in interpolation list:
      
   n <- node
   for ( i in 1:length(n$par) )
   {
      if ( ! is.na(n$names.par[i]) )
      {
         ind <- which(n$names.par[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter \"",
                    names(par)[ind[1]],"\"",sep="")
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$par[i] <- par[ind]
         }
      } 
   }
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.parfun1d <- function(x,
                                              attrib,   # data.frame, numeric
                                              par = NA,
                                              ...)
{
   node <- x
   
   # update parameters:

   n <- updatepar(node,par)
   
   # extract attributes:
   
   if ( is.data.frame(attrib) )
   {
      if ( length(which(names(attrib)==n$attrib)) != 1 )
      {
         warning("Node \"",node$name,"\": attribute \"",n$attrib,"\" not found",sep="")
         return(rep(NA,nrow(attrib)))
      }
      a <- attrib[,n$attrib]
   }
   else
   {
      if ( ! is.vector(attrib) )
      {
         warning("Node \"",node$name,"\": unknown format of attribute \"",n$attrib,"\"",sep="")
         return(NA)
      }
      if ( length(names(attrib)) == 0 )
      {
         a <- attrib
      }
      else
      {
         ind <- which(names(attrib)==n$attrib)
         if ( length(ind) != 1 )
         {
           if ( length(ind) > 1)
           {
              warning("Node \"",node$name,"\": multiple occurrences of attribute \"",
                      n$attrib,"\"",sep="")
           }
           else
           {
              warning("Node \"",node$name,"\": attribute \"",n$attrib,"\" not found",sep="")
           }
           return(NA)
         }
         a <- attrib[ind]
      }
   }
   
   # evaluate results:

   u <- do.call(n$name.fun,list(a,n$par))
   ind.out.of.range <- (a < n$range[1]) | (a > n$range[2])
   u <- ifelse(ind.out.of.range,NA,u)
   if ( sum(ind.out.of.range,na.rm=T) > 0 )
   {
      ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
      warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib,"\" out of range: ",
              paste(a[ind.not.na],collapse=","),sep="")
   }
   
   # return results:
   
   return(u)
}


# print:
# -----

print.utility.endnode.parfun1d <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.parfun1d <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   cat("attribute:      ",node$attrib,"\n")
   cat("attribute range:",node$range[1],"-",node$range[2],"\n")
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   cat("function type:  ",funtype,"\n")
   cat("required:       ",node$required,"\n")
   cat("function:       ",node$name.fun,"\n")
   cat("parameters:","\n")
   names.par <- ifelse(is.na(node$names.par),"",node$names.par)
   print(data.frame(names.par=names.par,par=node$par))
}


# plot:
# -----

plot.utility.endnode.parfun1d <- 
                   function(x,
                            par       = NA,
                            col       = utility.calc.colors(),
                            gridlines = c(0.2,0.4,0.6,0.8),
                            main      = "",
                            cex.main  = 1,
                            ...)
{
   node <- x
   n <- updatepar(node,par)
   utility.endnode.plot1d(node      = n,
                          col       = col,
                          gridlines = gridlines,
                          main      = main,
                          cex.main  = cex.main,
                          ...)
}


# ==============================================================================
# endnode for discrete factor attributes: 
# class "utility.endnode.discrete"
# ==============================================================================


# constructor:
# ------------

utility.endnode.discrete.create <- function(name.node,          # character(1)
                                            attrib.levels,      # data.frame
                                            u,                  # numeric(n)
                                            names.u     = rep(NA,length(u)),
                                            utility     = TRUE,
                                            required    = FALSE,
                                            col         = "black",
                                            shift.levels = 0)
{
   # consistency checks:

   check.ok <- T   
   if ( !is.data.frame(attrib.levels) )
   {
      cat("*** Warning: Attrib.levels must be a data frame","\n")
      check.ok <- F
   }
   if ( length(names(attrib.levels)) != length(unique(names(attrib.levels))) )
   {
      cat("*** Warning: Column names of attrib.levels must be different","\n")
      check.ok <- F
   }
   if ( nrow(attrib.levels) != length(u) )
   {
      cat("*** Warning: Number of rows of attrib.levels not equal to",
          "number of elements of u:",nrow(attrib.levels),length(u),"\n")
      check.ok <- F
   }
   if ( length(names.u) != length(u) )
   {
      cat("*** Warning: Number of elements of names.u not equal",
          "number of elements of u:",length(names.u),length(u),"\n")
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
   node$name          <- name.node
   node$description   <- "utility/value discrete attribute end node"
   node$type          <- "endnode"
   node$attrib.levels <- attrib.levels
   for ( i in 1:ncol(attrib.levels) ) 
   {
      node$attrib.levels[,i] <- as.character(node$attrib.levels[,i])
   }
   node$attrib        <- names(attrib.levels)
   node$u             <- u
   node$names.u       <- names.u
   node$required      <- required
   node$utility       <- utility
   node$col           <- col
   node$shift.levels   <- shift.levels
   class(node)        <- "utility.endnode.discrete" 
   
   # print and return class
   
   cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
   return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.discrete <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update adequate values in interpolation list:
      
   n <- node
   for ( i in 1:length(n$u) )
   {
      if ( ! is.na(n$names.u[i]) )
      {
         ind <- which(n$names.u[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter",
                    names(par)[ind[1]])
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$u[i] <- par[ind]
         }
      } 
   }
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.discrete <- function(x,
                                              attrib,   # data.frame
                                              par = NA,
                                              ...)
{
   node <- x
   
   # check availability of attributes:

   if ( ! is.data.frame(attrib) )
   {
      warning("Node \"",node$name,"\": attrib must be a data frame",sep="")
      return(NA)
   }
   ind <- match(node$attrib,names(attrib))
   if(sum(is.na(ind))>0)
   {
      ind.na <- is.na(ind)
      warning("Node \"",node$name,"\": attribute(s) \"",
              paste(node$attrib[ind.na],collapse=","),"\" not found",sep="")
      return(rep(NA,nrow(attrib)))
   }
   
   # check levels of attributes:
   
   for ( i in 1:ncol(node$attrib.levels) )
   {
      n <- names(node$attrib.levels)[i]                       # attribute name
      l <- unique(node$attrib.levels[,i]); l <- l[!is.na(l)]  # defined levels
      a <- unique(attrib[,n]); a <- a[!is.na(a)]              # requested levels
      ind.na <- is.na(match(a,l))
      if ( sum(ind.na) > 0 )
      {
         warning("Node \"",node$name,"\": unknown attribute level(s): \"",paste(a[ind.na],collapse=","),
                 "\" of attribute \"",n,"\"",sep="")
      }
   }
   
   # update parameters:

   n <- updatepar(node,par)
   
   # select rows compatible with conditioning attributes:

   u    <- rep(NA,nrow(attrib))
   calc <- rep(FALSE,nrow(attrib))
   for ( i in 1:ncol(n$attrib.levels) )  # evaluate NAs
   {
       calc <- calc | is.na(attrib[,names(n$attrib.levels)[i]])
   }
   while( TRUE )
   {
      # identify first row that has not yet been evaluated:
      
      startind <- match(FALSE,calc)
      
      # break if all were evaluated:
      
      if ( is.na(startind) ) break
      
      # find rows with the same attribute combinations:
      
      ind.attrib <- as.character(attrib[startind,names(n$attrib.levels)[1]]) == 
                    as.character(attrib[,names(n$attrib.levels)[1]])
      if ( ncol(n$attrib.levels) > 1 )
      {
         for ( i in 2:ncol(n$attrib.levels) )
         {
            ind.attrib <- 
               ind.attrib & 
               ( as.character(attrib[startind,names(n$attrib.levels)[i]]) == 
                 as.character(attrib[,names(n$attrib.levels)[i]]) )
         }
      }
      ind.attrib <- which(ind.attrib)
      
      # find corresponding value:
      
      ind.u <- as.character(attrib[startind,names(n$attrib.levels)[1]]) ==
               as.character(n$attrib.levels[,names(n$attrib.levels)[1]])
      if ( ncol(n$attrib.levels) > 1 )
      {
         for ( i in 2:ncol(n$attrib.levels) )
         {
            ind.u <- 
               ind.u & 
               ( as.character(n$attrib.levels[,names(n$attrib.levels)[i]]) == 
                 as.character(attrib[startind,names(n$attrib.levels)[i]]) )
         }
      }
      ind.u <- which(ind.u)
      
      # evaluate node for all attribute rows with same conditional values:
      
      if ( length(ind.u) == 1 )
      {
         u[ind.attrib] <- n$u[ind.u]
      }
      else
      {
         if ( length(ind.u) > 1 )
         {
            warning("Node \"",node$name,"\": multiple combinations of the same",
                    "attribute levels in node \"",n$name,"\"",sep="")
         }
      }            
      calc[ind.attrib] <- T
   }      
   
   # return results:
   
   return(u)
}


# print:
# -----

print.utility.endnode.discrete <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.discrete <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   cat("attribute(s):   ",paste(node$attrib,collapse=","),"\n")
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   cat("function type:  ",funtype,"\n")
   cat("required:       ",node$required,"\n")
   cat("attribute/value combinations:","\n")
   names.u <- ifelse(is.na(node$names.u),"",node$names.u)
   print(cbind(node$attrib.levels,u=node$u,names.u=names.u))
}


# plot:
# -----

plot.utility.endnode.discrete <- 
                         function(x,
                                  par       = NA,
                                  col       = utility.calc.colors(),
                                  gridlines = c(0.2,0.4,0.6,0.8),
                                  main      = "",
                                  cex.main  = 1,
                                  ...)
{
   node <- x
   length = 101
   n <- updatepar(node,par)
   title <- main; if ( nchar(title) == 0 ) title <- n$name
   funtype <- "utility"; if ( !n$utility ) funtype <- "value"
   plot(numeric(0),numeric(0),type="l",
        xlim=c(0,1),ylim=c(0,1),
        xlab=paste(n$attrib,collapse=","),ylab=funtype,main=title,
        xaxs="i",yaxs="i",xaxt="n",cex.main=cex.main,...)
   labels=character(length(n$u))
   for ( i in 1:length(n$u) )
   {
      labels[i] <- paste(as.character(n$attrib.levels[i,]),collapse=",")
   }
   axis(side=1,at=((1:length(n$u))-0.5)/length(n$u),labels=labels)
   
   # color axis and plot gridlines:
   
   if ( length(col)>1 & !node$utility )
   {
      du <- 1/(length-1)
      midpoints <- seq(du,1-du,length=length-1)
      cols <- utility.get.colors(midpoints,col)
      for ( i in 1:(length-1) )
      {
         lines(c(1,1)*0.001,
               c((i-1)*du,i*du),col=cols[i],lwd=3)
      }
      if ( ! is.na(gridlines[1]) )
      {
         for ( level in gridlines ) abline(h=level,lty="dashed")
      }
   }
   
   # plot points:
   
   color <- "black"; if(length(col)==1) color <- col
   points(((1:length(n$u))-0.5)/length(n$u),n$u,pch=19,col=color,xpd=TRUE)
}


# ==============================================================================
# endnode for combining other endnodes conditional on factor attributes: 
# class "utility.endnode.cond"
# ==============================================================================


# constructor:
# ------------

utility.endnode.cond.create <- function(name.node,          # character(1)
                                        attrib.levels,      # data.frame
                                        nodes,              # list of nodes
                                        utility     = TRUE,
                                        required    = FALSE,
                                        col         = "black",
                                        shift.levels = 0)
{
   # consistency checks:

   check.ok <- T   
   if ( !is.data.frame(attrib.levels) )
   {
      cat("*** Warning: attrib.levels must be a data frame","\n")
      check.ok <- F
   }
   if ( length(names(attrib.levels)) != length(unique(names(attrib.levels))) )
   {
      cat("*** Warning: cColumn names of attrib.levels must be different","\n")
      check.ok <- F
   }
   if ( nrow(attrib.levels) != length(nodes) )
   {
      cat("*** Warning: Number of rows of attrib.levels not equal to",
          "number of nodes provided:",nrow(attrib.levels),length(nodes),"\n")
      check.ok <- F
   }
   if ( length(nodes) < 1 )
   {
      cat("*** Warning: No nodes provided","\n")
      check.ok <- F
   }
   for ( i in 1:length(nodes) )
   {
      if ( nodes[[i]]$utility != utility )
      {
         funtype   <- "utility"; if ( !utility )            funtype   <- "value"
         funtype.i <- "utility"; if ( !nodes[[i]]$utility ) funtype.i <- "value"
         cat("***Warning: incompatible function types: new node is of type",
             funtype,"node",nodes[[i]]$name," is of type",funtype.i,"\n")
         check.ok <- F
      }
   }
   if ( ! check.ok )
   {
      cat("*** Warning: Node \"",name.node,"\" could not be constructed","\n",
          sep="")
      return(NA)
   }

   # construct class:
   
   node <- list()
   node$name          <- name.node
   node$description   <- "utility/value conditional combination end node"
   node$type          <- "endnode"
   node$attrib.levels <- attrib.levels
   for ( i in 1:ncol(attrib.levels) ) 
   {
      node$attrib.levels[,i] <- as.character(node$attrib.levels[,i])
   }
   node$attrib        <- names(attrib.levels)
   for ( i in 1:length(nodes) )
   {
      node$attrib     <- c(node$attrib,nodes[[i]]$attrib)
   }
   node$attrib        <- unique(node$attrib)
   node$nodes         <- nodes
   node$required      <- required
   node$utility       <- utility
   node$col           <- col
   node$shift.levels   <- shift.levels
   class(node)        <- "utility.endnode.cond" 
   
   # print return class
   
   cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
   return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.cond <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update conditional nodes:

   n <- node
   for ( i in 1:length(n$nodes) )
   {
      n$nodes[[i]] <- updatepar(n$nodes[[i]],par)
   }
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.cond <- function(x,
                                          attrib,   # data.frame
                                          par=NA,
                                          ...)
{
   node <- x
   
   # check availability of attributes:

   if ( ! is.data.frame(attrib) )
   {
      warning("Node \"",node$name,"\": attrib must be a data frame",sep="")
      return(NA)
   }
   ind <- match(node$attrib,names(attrib))
   if(sum(is.na(ind))>0)
   {
      ind.na <- is.na(ind)
      warning("Node \"",node$name,"\": attribute(s) \"",
              paste(node$attrib[ind.na],collapse=","),"\" not found",sep="")
      return(rep(NA,nrow(attrib)))
   }

   # update parameters:

   n <- updatepar(node,par)
   
   # select rows compatible with conditioning attributes:
   
   u    <- rep(NA,nrow(attrib))
   calc <- rep(FALSE,nrow(attrib))
   for ( i in 1:ncol(n$attrib.levels) )  # evaluate NAs
   {
       calc <- calc | is.na(attrib[,names(n$attrib.levels)[i]])
   }
   while( TRUE )
   {
      # identify first row that has not yet been evaluated:
      
      startind <- match(FALSE,calc)
      
      # break if all were evaluated:
      
      if ( is.na(startind) ) break
      
      # find rows with the same attribute combinations:
      
      ind.attrib <- as.character(attrib[startind,names(n$attrib.levels)[1]]) == 
                    as.character(attrib[,names(n$attrib.levels)[1]])
      if ( ncol(n$attrib.levels) > 1 )
      {
         for ( i in 2:ncol(n$attrib.levels) )
         {
            ind.attrib <- 
               ind.attrib & 
               ( as.character(attrib[startind,names(n$attrib.levels)[i]]) == 
                 as.character(attrib[,names(n$attrib.levels)[i]]) )
         }
      }
      ind.attrib <- which(ind.attrib)
      
      # find corresponding node:

      ind.node <- as.character(attrib[startind,names(n$attrib.levels)[1]]) ==
                  as.character(n$attrib.levels[,names(n$attrib.levels)[1]])
      if ( ncol(n$attrib.levels) > 1 )
      {
         for ( i in 2:ncol(n$attrib.levels) )
         {
            ind.node <- 
               ind.node & 
               ( as.character(n$attrib.levels[,names(n$attrib.levels)[i]]) == 
                 as.character(attrib[startind,names(n$attrib.levels)[i]]) )
         }
      }
      ind.node <- which(ind.node)
      
      # evaluate node for all attribute rows with same conditional values:

      if ( length(ind.node) > 0 )
      {      
         u[ind.attrib] <- evaluate(n$nodes[[ind.node[1]]],attrib[ind.attrib,])
         if ( length(ind.node) > 1 )
         {
            cat("*** Warning: multiple combinations of the same",
                "attribute levels in node",n$name,"\n")
         }
      }            
      calc[ind.attrib] <- T
   }      
   
   # return results:
   
   return(u)
}


# print:
# -----

print.utility.endnode.cond <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.cond <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   cat("function type:  ",funtype,"\n")
   cat("required:       ",node$required,"\n")
   cat("attribute/node combinations:","\n")
   nodes.names <- character(0)
   for ( i in 1:length(node$nodes) ) nodes.names[i] <- node$nodes[[i]]$name
   print(cbind(node$attrib.levels,node=nodes.names))
   for ( i in 1:length(node$nodes) ) 
   {
      cat("**","\n")
      summary(node$nodes[[i]])
   }
}


# plot:
# -----

plot.utility.endnode.cond <-
                         function(x,
                                  par       = NA,
                                  col       = utility.calc.colors(),
                                  gridlines = c(0.2,0.4,0.6,0.8),
                                  main      = "",
                                  cex.main  = 1,
                                  nodes     = x$name,
                                  ...)
{
   node <- x
   if ( is.na(nodes[1]) | ! is.na(match(node$name,nodes)) )
   {
      nrow <- floor(sqrt(length(node$nodes)))
      ncol <- floor(length(node$nodes)/nrow+0.999)
      par.def <- par(no.readonly=T)
      par(mfrow=c(nrow,ncol),mar=c(4.3,3.8,2.8,0.8),oma=c(0,0,2,0)) 
      for ( i in 1:length(node$nodes) )             # c(bottom, left, top, right)
      {
         title <- main
         for ( j in 1:ncol(node$attrib.levels) )
         {
            title <- paste(title," ",colnames(node$attrib.levels)[j],"=",
                           as.character(node$attrib.levels[i,j]),sep="")
         }
         plot(node$nodes[[i]],par=par,col=col,gridlines=gridlines,main=title,cex.main=cex.main,...)
      }
      mtext(node$name,outer=TRUE,cex=cex.main)
      par(par.def)
   }
   if ( length(node$nodes) > 0 )
   {
      for ( i in 1:length(node$nodes) )
      {
         if ( is.na(nodes[1]) | !is.na(match(node$nodes[[i]]$name,nodes)) )
         {
            plot(node$nodes[[i]],
                 par=par,
                 col=col,
                 gridlines=gridlines,
                 cex.main=cex.main,
                 ...)
         }
      }
   }
}


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
                           shift.levels  = 0)
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
   class(node)       <- "utility.aggregation" 
   
   # return class
   
   cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
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
   u.agg <- apply(u.agg.input,1,n$name.fun,n$par)
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
                ...)
}


# ==============================================================================
# conversion node from values to utilities with interpolation: 
# class "utility.conversion.intpol"
# ==============================================================================


# constructor:
# ------------

utility.conversion.intpol.create <- function(name.node,    # character(1)
                                             node,         # character(1)
                                             x,            # numeric(n)
                                             u,            # numeric(n)
                                             names.x     = rep(NA,length(x)),
                                             names.u     = rep(NA,length(u)),
                                             required    = FALSE,
                                             col         = "black",
                                             shift.levels = 0)
{
   # consistency checks:

   check.ok <- T   
   if ( length(x) != length(u) )
   {
      cat("*** Warning: x and u of different length:",
          length(x),length(u),"\n")
      check.ok <- F
   }
   if ( length(names.x) != length(names.u) )
   {
      cat("*** Warning: names.x and names.u of different length:",
          length(names.x),length(names.u),"\n")
      check.ok <- F
   }
   if ( length(x) != length(names.x) )
   {
      cat("*** Warning: x and names.x of different length:",
          length(x),length(names.x),"\n")
      check.ok <- F
   }
   if ( ! utility.check.name(name.node,node) )
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
   
   n <- list()
   n$name         <- name.node
   n$description  <- "utility/value interpolation conversion node"
   n$type         <- "conversionnode"
   n$nodes        <- list(node)
   n$x            <- x
   n$u            <- u
   n$names.x      <- names.x
   n$names.u      <- names.u
   n$required     <- required
   n$num.required <- 1
   n$utility      <- TRUE
   n$col          <- col
   n$shift.levels  <- shift.levels
   class(n)       <- "utility.conversion.intpol" 
   
   # print and return class

   cat(n$description," \"",name.node,"\" constructed","\n",sep="")   
   return(n)
}


# update parameter values:
# ------------------------

updatepar.utility.conversion.intpol <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update adequate values in interpolation list:
      
   n <- node
   for ( i in 1:length(n$x) )
   {
      if ( ! is.na(n$names.x[i]) )
      {
         ind <- which(n$names.x[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter",
                    names(par)[ind[1]],sep="")
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$x[i] <- par[ind]
         }
      } 
      if ( ! is.na(n$names.u[i]) )
      {
         ind <- which(n$names.u[i] == names(par) )
         if ( length(ind) > 1 )
         {
            warning("Node \"",node$name,"\": multiple occurrences of parameter",
                    names(par)[ind[1]],sep="")
            ind <- ind[1]
         }
         if ( length(ind) == 1 )
         {
            n$u[i] <- par[ind]
         }
      } 
   }
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.cond.utility.conversion.intpol <- function(x,v,...)
{
   node <- x
   u <- approx(x=node$x,y=node$u,xout=v)$y
   return(u)
}


evaluate.utility.conversion.intpol <- function(x,
                                               attrib,   # data.frame, numeric
                                               par = NA,
                                               ...)
{
   node <- x
   
   # update parameters:

   n <- updatepar(node,par)
   
   # evaluate results:
   
   v <- evaluate(n$nodes[[1]],attrib)
   if ( ! is.data.frame(v) )
   {
      v <- as.data.frame(v)
   }
   u <- evaluate.cond(n,v[,1])
   ind <- !is.na(u) & (u<0 | u>1)
   if ( sum(ind) > 0 )
   {
      warning("Node \"",node$name,"\": node \"",n$name,"\" produced values outside of [0,1]: ",
              paste(u[ind],collapse=","),sep="")
   }
   u <- as.data.frame(u)
   names(u) <- node$name
   
   # return results:
   
   u <- cbind(u,v)
   rownames(u) <- rownames(attrib)

   return(u)
}


# print:
# -----

print.utility.conversion.intpol <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.conversion.intpol <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   cat("function type:  ","utility","\n")
   cat("required:       ",node$required,"\n")
   cat("data pairs:","\n")
   names.x <- ifelse(is.na(node$names.x),"",node$names.x)
   names.u <- ifelse(is.na(node$names.u),"",node$names.u)
   print(data.frame(names.x=names.x,x=node$x,u=node$u,names.u=names.u))
   for ( i in 1:length(node$nodes) ) 
   {
      cat("***","\n")
      summary(node$nodes[[i]])
   }
}


# plot:
# -----

plot.utility.conversion.intpol <- 
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
                ...)
}


# ==============================================================================
# conversion node from values to utilities with parametric function: 
# class "utility.conversion.parfun"
# ==============================================================================


# constructor:
# ------------

utility.conversion.parfun.create <- function(name.node,    # character(1)
                                             node,         # node
                                             name.fun,     # name of f(a,par)
                                             par,          # numeric(n)
                                             names.par   = rep(NA,length(par)),
                                             required    = FALSE,
                                             col         = "black",
                                             shift.levels = 0)
{
   # consistency checks:

   check.ok <- T   
   if ( length(par) != length(names.par) )
   {
      cat("*** Warning: par and names.par of different length:",
          length(par),length(names.par),"\n")
      check.ok <- F
   }
   if ( ! utility.check.name(name.node,list(node)) )
   {
      cat("*** Warning: node with same name \"",name.node,"\" exists already ",
          "as sub-node","\n")
      check.ok <- F
   }
   if ( ! check.ok )
   {
      cat("*** Warning: node \"",name.node,"\" could not be constructed","\n",
          sep="")
      return(NA)
   }

   # construct class:
   
   n <- list()
   n$name         <- name.node
   n$description  <- "utility/value parametric function conversion node"
   n$type         <- "utility.conversion.parfun"
   n$nodes        <- list(node)
   n$name.fun     <- name.fun
   n$par          <- par
   n$names.par    <- names.par
   n$required     <- required
   n$num.required <- 1
   n$utility      <- TRUE
   n$col          <- col
   n$shift.levels  <- shift.levels
   class(n)       <- "utility.conversion.parfun" 
   
   # print and return class
   
   cat(n$description," \"",name.node,"\" constructed","\n",sep="")   
   return(n)
}


# update parameter values:
# ------------------------

updatepar.utility.conversion.parfun <- function(x,par=NA,...)
{
   node <- x
   
   # check availabiliy of named parameter vector:
   
   if ( length(names(par)) == 0 ) return(node)

   # update adequate values in interpolation list:
      
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
   
   # return updated node:
   
   return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.cond.utility.conversion.parfun <- function(x,v,...)
{
   node <- x
   u <- do.call(node$name.fun,list(v,node$par))
   return(u)
}


evaluate.utility.conversion.parfun <- function(x,
                                               attrib,   # data.frame, numeric
                                               par = NA,
                                               ...)
{
   node <- x
   
   # update parameters:

   n <- updatepar(node,par)
   
   # evaluate results:
   
   v <- evaluate(n$nodes[[1]],attrib)
   if ( ! is.data.frame(v) )
   {
      v <- as.data.frame(v)
   }
   u <- evaluate.cond(n,v[,1])
   u <- as.data.frame(u)
   names(u) <- n$name
   ind <- !is.na(u) & (u<0 | u>1)
   if ( sum(ind) > 0 )
   {
     warning("Node \"",node$name,"\": node \"",n$name,"\" produced values outside of [0,1]: ",
             paste(u[ind],collapse=","),sep="")
   }
   
   # return results:
   
   u <- cbind(u,v)
   rownames(u) <- rownames(attrib)
   
   # return results:
   
   return(u)
}


# print:
# -----

print.utility.conversion.parfun <- function(x,...)
{
   cat(paste(rep("-",50),collapse=""),"\n")
   summary(x,...)
   cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.conversion.parfun <- function(object,...)
{
   node <- object
   cat(node$name,"\n")
   cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
   cat(node$description,"\n")
   cat("node     :      ",node$nodes[[1]]$name,"\n")
   cat("function type:  ","utility","\n")
   cat("required:       ",node$required,"\n")
   cat("function:       ",node$name.fun,"\n")
   cat("parameters:","\n")
   names.par <- ifelse(is.na(node$names.par),"",node$names.par)
   print(data.frame(names.par=names.par,par=node$par))
   for ( i in 1:length(node$nodes) ) 
   {
      cat("***","\n")
      summary(node$nodes[[i]])
   }
}


# plot:
# -----

plot.utility.conversion.parfun <- 
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
                cex.main    = cex.main,
                cex.nodes   = cex.nodes,
                cex.attrib  = cex.attrib,
                f.reaches   = f.reaches,
                f.nodes     = f.nodes,
                with.attrib = with.attrib,
                ...)
}


# ==============================================================================





