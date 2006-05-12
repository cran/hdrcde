#modalreg 01/05/2006 JE
# modified 11/5/06 RH

modalreg <- function(x, y, xfix=seq(min(x),max(x),l=50),  a, b, deg = 0, iter= 30, P = 2,
  start="e", shrink=TRUE, shrink.const=10, plot.type=c("p",1),
  labels= c("","x","y"), pch=".", ...)
{
    require(hdrcde)

    if (missing(xfix))
        xfix <- x
    if (missing(a) || missing(b))
    {
        h.nr <- cde.bandwidths(x, y,method=1, deg=deg, ...)
        a <- h.nr$a
        b <- h.nr$b
    }

    if (P ==1)
        ynull <- quantile(y,probs=0.5)
    else
    {
         if (start=="q")
            ynull <- quantile(y, probs=seq(0,1,by=1/(P-1)))
         else if (start=="e" || start=="v")
            ynull <- seq(min(y),max(y),length=P)
         else
            ynull <- runif(P,min(y),max(y))
    }

    n <- length(x)
    save.regression <- matrix(0, P, length(xfix))

    Alphabet <- c("A","B", "C", "D", "E", "F","G","H")
    alphabet <- c("+","x", "*", "d", "e", "f","g","h")

    if (plot.type[1]!="n")
    {
        plot(x,y,pch=pch,main=labels[1],xlab=labels[2],ylab=labels[3], col='grey')
        if (plot.type[2]==1)
            points(rep(min(x),P),ynull,col=2:(P+1),pch=Alphabet[1:P])
    }


    span.area <- (max(x)-min(x))*(max(y)-min(y))

    for (i in 1: length(xfix))
    {
        for (p in 1:P)
        {
            if (start!="v" || i==1)
                current.regression <- ynull[p]
            else
                current.regression <-  save.regression[p,i-1]
            old.regression <- -1000
            for (j in 1:iter)
            {
                old.regression <- current.regression
                if (deg==1)
                    current.regression <- cond.linear.meanshift(x,y, xfix[i],current.regression,a,b)
                else if (deg==0)
                    current.regression <- cond.meanshift(x,y, xfix[i],current.regression,a,b)
                if (current.regression=="NaN")
                {
                    current.regression <- old.regression
                    break()
                }
            }
            save.regression[p,i] <- current.regression
        }
        if (i%%10 ==0)
        cat(i,".." )
    }
    cat("\n")

    kde <- matrix(0,P,length(xfix))
    if (shrink==TRUE)
    {
        for (i in 1: length(xfix))
        {
            for (p in 1:P)
            {
                escape <- ifelse(p==1,2,p-1)
                kde[p,i] <-  kde2d.point(x,y,xfix[i],save.regression[p,i],a,b)
                if (kde[p,i]<1/(shrink.const*span.area))
                    save.regression[p,i] <- save.regression[escape,i]
            }
        }
    }

    if (plot.type[1]!="n")
    {
        for (p in 1:P)
        {
            if (plot.type[1]=="p")
                points(xfix, save.regression[p,],cex=1,col=p+1,pch=alphabet[p])
            else
                lines(xfix, save.regression[p,],cex=2)
        }
    }
    h <- c(a,b); names(h) <- c("a","b")
    list("Fitted.values"=save.regression, "Bandwidths"=h, "Density"= kde, "Threshold" = 1/(shrink.const*span.area))
}


######### Auxiliary functions

cond.meanshift <- function(x,y,x0, y0, a, b)
{
    sum(kern(x,x0,a)*gern(y,y0,b)*y)/(sum(kern(x,x0,a)*gern(y,y0,b)))
}

cond.linear.meanshift <- function(x,y,x0, y0, a, b)
{
    sn1 <- sum(kern(x,x0,a)*(x0-x))
    sn2 <- sum(kern(x,x0,a)*(x0-x)^2)
    sum(kern(x,x0,a)*(sn2-(x0-x)*sn1)*gern(y,y0,b)*y)/(sum(kern(x,x0,a)*(sn2-(x0-x)*sn1)*gern(y,y0,b)))
}

#Kernel function K1 (horizontal, "kern")
kern <-  function(x, x0 = 0, h = 1)
{
      1/h * dnorm((x0 - x)/h)
}

#Kernel function G (vertical, "gern"), is equal to K2 if G Gaussian.
gern <- kern

#Fast kernel density estimate   (faster than kde in package ks)
kde2d.point <- function(x, y, x0, y0, a, b)
{
    1/(length(x))*sum(kern(x,x0,a)*gern(y,y0,b))
}


