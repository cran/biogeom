\name{Neocinnamomum}
\docType{data}
\alias{Neocinnamomum}
\title{
Leaf Boundary Data of Seven Species of \emph{Neocinnamomum}
}
\description{
The data consist of the leaf boundary data of seven species of \emph{Neocinnamomum}.
}
\usage{
data(Neocinnamomum)
}
\details{
In the data set, there are four columns of variables: \code{Code}, \code{LatinName}, \code{x}, and \code{y}. 
\code{Code} saves the codes of individual leaves; 
\code{LatinName} saves the Latin names of the seven species of \emph{Neocinnamomum}; 
\code{x} saves the \eqn{x} coordinates of the leaf boundary in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the leaf boundary in the Cartesian coordinate system (cm).
}
\references{
Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Yu, K., Niklas, K.J., Schrader, J., Song, Y., Zhu, R., Li, Y., Wei, H., Ratkowsky, D.A. (2021)  
A general model for describing the ovate leaf shape. \emph{Symmetry}, 13, 1524. \doi{10.3390/sym13081524}
}
\examples{
data(Neocinnamomum)

uni.C <- sort( unique(Neocinnamomum$Code) )
ind   <- 2
Data  <- Neocinnamomum[Neocinnamomum$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y
length(x0)

Res1  <- adjdata(x0, y0, ub.np=200, len.pro=1/20)
x1    <- Res1$x
y1    <- Res1$y
length(x1)

dev.new()
plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) ) 

}
