\name{kp}
\docType{data}
\alias{kp}
\title{
Boundary Data of the Vertical Projections of \emph{Koelreuteria paniculata} Fruit
}
\description{
The data consist of the boundary data of four vertical projections of \emph{K. paniculata} fruit  
sampled at Nanjing Forestry University campus in early October 2021.
}
\usage{
data(kp)
}
\details{
In the data set, there are three columns of variables: \code{Code}, \code{x}, and \code{y}. 
\code{Code} saves the codes of individual fruit; 
\code{x} saves the \eqn{x} coordinates of the vertical projections of fruit in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the vertical projections of fruit in the Cartesian coordinate system (cm).
}
\references{
Li, Y., Quinn, B.K., Gielis, J., Li, Y., Shi, P. (2022) Evidence that supertriangles exist 
in nature from the vertical projections of \emph{Koelreuteria paniculata} fruit. \emph{Symmetry} 
14, 23. \doi{10.3390/sym14010023}
}
\examples{
data(kp)

uni.C <- sort( unique(kp$Code) )
ind   <- 1
Data  <- kp[kp$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1 <- adjdata(x0, y0, ub.np=200, len.pro=1/20)
dev.new()
plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

\donttest{
  x1        <- Res1$x
  y1        <- Res1$y
  x0.ini    <- mean( x1 )
  y0.ini    <- mean( y1 )
  theta.ini <- pi
  a.ini     <- 0.9
  n1.ini    <- c(1, 4)
  n2.ini    <- 5
  n3.ini    <- c(5, 10, 15)
  ini.val   <- list(x0.ini, y0.ini, theta.ini, 
                    a.ini, n1.ini, n2.ini, n3.ini)

  Res2 <- fitGE( GE, x=x1, y=y1, ini.val=ini.val, 
                 m=3, simpver=5, nval=1, unit="cm",  
                 par.list=FALSE, fig.opt=TRUE, angle=NULL, 
                 control=list(reltol=1e-20, maxit=20000), 
                 np=2000 )
}

}
