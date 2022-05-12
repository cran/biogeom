\name{starfish}
\docType{data}
\alias{starfish}
\title{
Boundary Data of Eight Sea Stars
}
\description{
The data consist of the boundary data of eight sea stars from five species.
}
\usage{
data(starfish)
}
\details{
In the data set, there are four columns of variables: \code{Code}, \code{LatinName}, \code{x}, and \code{y}. 
\code{Code} saves the codes of individual sea stars; 
\code{LatinName} saves the Latin names of the eight sea stars; 
\code{x} saves the \eqn{x} coordinates of the eight sea stars in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the eight sea stars in the Cartesian coordinate system (cm). 
In \code{Code}, codes 1-9 represent \emph{Anthenoides tenuis}, \emph{Culcita schmideliana} sample 1, 
\emph{Culcita schmideliana} sample 2, \emph{Culcita schmideliana} sample 3, \emph{Stellaster equestris}, \emph{Tosia australis}, 
\emph{Tosia magnifica} sample 1, and \emph{Tosia magnifica} sample 2, respectively. See Table A1 published in Shi et al. (2020).
}
\references{
Shi, P., Ratkowsky, D.A., Gielis, J. (2020) The generalized Gielis geometric 
equation and its application. \emph{Symmetry} 12, 645. \doi{10.3390/sym12040645}
}
\examples{
data(starfish)

uni.C <- sort( unique(starfish$Code) )
ind   <- 2
Data  <- starfish[starfish$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

dev.new()
plot( x0, y0, asp=1, type="l", cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )
length(x0)

Res1 <- adjdata(x0, y0, ub.np=400, times=1.2, len.pro=1/20)
x1   <- Res1$x
y1   <- Res1$y

dev.new()
plot( x1, y1, asp=1, type="l", cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )
length(x1)

}
