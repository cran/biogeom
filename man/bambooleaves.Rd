\name{bambooleaves}
\docType{data}
\alias{bambooleaves}
\title{
Leaf Boundary Data of \emph{Phyllostachys incarnata} T. H. Wen (Poaceae: Bambusoideae)
}
\description{
The data consist of the boundary data of six leaves of \emph{P. incarnata} 
sampled at Nanjing Forestry University campus in early December 2016.
}
\usage{
data(bambooleaves)
}
\details{
In the data set, there are three columns of variables: \code{Code}, \code{x}, and \code{y}. 
\code{Code} saves the codes of individual leaves; 
\code{x} saves the \eqn{x} coordinates of the leaf boundary in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the leaf boundary in the Cartesian coordinate system (cm).
}
\references{
Lin, S., Shao, L., Hui, C., Song, Y., Reddy, G.V.P., Gielis, J., Li, F., Ding, Y., Wei, Q., 
Shi, P. (2018) Why does not the leaf weight-area allometry of bamboos follow the 3/2-power law? 
\emph{Frontiers in Plant Science} 9, 583. \doi{10.3389/fpls.2018.00583}

Shi, P., Ratkowsky, D.A., Li, Y., Zhang, L., Lin, S., Gielis, J. (2018) 
General leaf-area geometric formula exists for plants - Evidence from the simplified Gielis equation. 
\emph{Forests} 9, 714. \doi{10.3390/f9110714}
}
\examples{
data(bambooleaves)

uni.C <- sort( unique(bambooleaves$Code) )
ind   <- 1
Data  <- bambooleaves[bambooleaves$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

dev.new()
plot( x0, y0, asp=1, type="l", cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )
length(x0)

Res1 <- adjdata(x0, y0, ub.np=600, len.pro=1/20)
dev.new()
plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

graphics.off()
}
