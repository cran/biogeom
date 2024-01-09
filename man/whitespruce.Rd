\name{whitespruce}
\encoding{UTF-8}
\docType{data}
\alias{whitespruce}
\title{
Planar Coordinates of \emph{Picea glauca} Tree Rings
}
\description{
The data consist of the planar coordinates of \emph{Picea glauca} tree rings.
}
\usage{
data(whitespruce)
}
\details{
In the data set, there are three columns of variables: \code{Code}, \code{x}, and \code{y}. 
\code{Code} saves the age codes of tree rings from the 2nd year to the 44th year; 
\code{x} saves the \eqn{x} coordinates of the tree rings in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the tree rings in the Cartesian coordinate system (cm).
}
\references{
Shi, P., Huang, J., Hui, C., Grissino-Mayer, H.D., Tardif, J., Zhai, L., Wang, F., Li, B. (2015) 
Capturing spiral radial growth of conifers using the superellipse to model tree-ring geometric shape. 
\emph{Frontiers in Plant Science} 6, 856. \doi{10.3389/fpls.2015.00856}
}
\examples{
data(whitespruce)

uni.C <- sort( unique(whitespruce$Code) )  
Data  <- whitespruce[whitespruce$Code==uni.C[10], ]
x0    <- Data$x
y0    <- Data$y
Res1  <- adjdata(x0, y0, ub.np=2000, len.pro=1/20)

plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
      xlim=c(3, 13), ylim=c(3, 13), col="grey73", lwd=2,
      xlab=expression(italic(x)), ylab=expression(italic(y)) )

\donttest{

  uni.C <- sort( unique(whitespruce$Code) )  
  for(i in 1:length(uni.C)){
    Data  <- whitespruce[whitespruce$Code==uni.C[i], ]
    x0    <- Data$x
    y0    <- Data$y

    Res1  <- adjdata(x0, y0, ub.np=200, len.pro=1/10)
  
    if(i == 1){
      plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
            xlim=c(3, 13), ylim=c(3, 13), col=1, lwd=1,
            xlab=expression(italic(x)), ylab=expression(italic(y)) )
    }
    if(i > 1) lines(Res1$x, Res1$y, col=1, lwd=1)
   
  }

  uni.C   <- sort( unique(whitespruce$Code) ) 
  uni.C   <- uni.C[1:12]
  Length  <- c()
  results <- data.frame(Code=c(), x0=c(), y0=c(), theta=c(),
                   a=c(), k=c(), n1=c(), r.sq=c(), RSS=c(), N=c())
  for(i in 1:length(uni.C)){
    Data      <- whitespruce[whitespruce$Code==uni.C[i], ]
    x0        <- Data$x
    y0        <- Data$y
    Res1      <- adjdata(x0, y0, ub.np=200, len.pro=1/10)
    x1        <- Res1$x
    y1        <- Res1$y
    x0.ini    <- mean( x1 )
    y0.ini    <- mean( y1 )
    theta.ini <- c(0, pi/4, pi/2)
    a.ini     <- 0.9
    k.ini     <- 1
    n1.ini    <- c(1.5, 2, 2.5)
    ini.val   <- list(x0.ini, y0.ini, theta.ini, 
                      a.ini, k.ini, n1.ini)

    print(paste("Progress: ", i, "/", length(uni.C), sep=""))
    H <- NULL
    try( H <- fitGE(GE, x=x1, y=y1, ini.val=ini.val, 
                    m=4, simpver=9, unit="cm", par.list=FALSE,
                    stand.fig=FALSE, angle=NULL, fig.opt=FALSE,  
                    control=list(reltol=1e-20, maxit=20000), 
                    np=2000), silent=TRUE )
    if(is.null(H)){
      RE <- data.frame(Code=uni.C[i], x0=NA, y0=NA, theta=NA,
               a=NA, k=NA, n1=NA, r.sq=NA, RSS=NA, N=NA)
    }
    if(!is.null(H)){
      RE     <- data.frame(Code=uni.C[i], x0=H$par[1], y0=H$par[2], 
                  theta=H$par[3], a=H$par[4], k=H$par[5], n1=H$par[6],
                  r.sq=H$r.sq, RSS=H$RSS, N=H$sample.size)
      Length <- c(Length, max(max(H$y.stand.pred)[1]-min(H$y.stand.pred)[1], 
                  max(H$x.stand.pred)[1]-min(H$x.stand.pred)[1])[1])
      if(i == 1){
        plot(H$x.obs, H$y.obs, asp=1, xlim=c(7.4, 8.6), ylim=c(7.4, 8.6),
             cex.lab=1.5, cex.axis=1.5, type="l", lwd=2, col="grey70",
             xlab=expression(italic(x)), ylab=expression(italic(y)))
        lines(H$x.pred, H$y.pred, col=2)         
      } 
      if(i > 1){
        lines(H$x.obs, H$y.obs, lwd=2, col="grey70")
        lines(H$x.pred, H$y.pred, col=2)   
      }
      
    }
    results <- rbind(results, RE)
  }

  # To adjust the estimates of partial parameters to ensure k <= 1
  results2            <- results
  Ind                 <- results$k > 1
  results2$theta[Ind] <- results$theta[Ind] + pi/2
  results2$a[Ind]     <- results$a[Ind] * results$k[Ind]^(1/results$n1[Ind])
  results2$k[Ind]     <- 1/results$k[Ind] 
  results2
  Length/2

}  

}
