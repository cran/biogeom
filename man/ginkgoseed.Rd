\name{ginkgoseed}
\docType{data}
\alias{ginkgoseed}
\title{
Boundary Data of the Side Projections of \emph{Ginkgo biloba} Seeds
}
\description{
The data consist of the boundary data of four side projections of \emph{G. biloba} (Cultivar 'Fozhi') seeds 
sampled at Nanjing Forestry University campus on September 23, 2021.
}
\usage{
data(ginkgoseed)
}
\details{
In the data set, there are three columns of variables: \code{Code}, \code{x}, and \code{y}. 
\code{Code} saves the codes of individual fruit; 
\code{x} saves the \eqn{x} coordinates of the side projections of seeds in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the side projections of seeds in the Cartesian coordinate system (cm).
}
\references{
Tian, F., Wang, Y., Sandhu, H.S., Gielis, J., Shi, P. (2020) Comparison of seed morphology of two ginkgo cultivars. 
\emph{Journal of Forestry Research} 31, 751\eqn{-}758. \doi{10.1007/s11676-018-0770-y}
}
\examples{
data(ginkgoseed)

uni.C <- sort( unique(ginkgoseed$Code) )
ind   <- 1
Data  <- ginkgoseed[ginkgoseed$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1 <- adjdata(x0, y0, ub.np=2000, len.pro=1/20)
dev.new()
plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

\donttest{
  x1        <- Res1$x
  y1        <- Res1$y
  x0.ini    <- mean( x1 )
  y0.ini    <- mean( y1 )
  theta.ini <- pi/4
  a.ini     <- 1
  n1.ini    <- seq(0.6, 1, by=0.1)
  n2.ini    <- 1
  n3.ini    <- 1
  ini.val   <- list(x0.ini, y0.ini, theta.ini, 
                    a.ini, n1.ini, n2.ini, n3.ini)

  Res2 <- fitGE( GE, x=x1, y=y1, ini.val=ini.val, 
                 m=2, simpver=5, nval=1, unit="cm",  
                 par.list=FALSE, fig.opt=TRUE, angle=NULL, 
                 control=list(reltol=1e-20, maxit=20000), 
                 np=2000 )
}

}
