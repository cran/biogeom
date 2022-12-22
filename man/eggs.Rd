\name{eggs}
\docType{data}
\alias{eggs}
\title{
Egg Boundary Data of Nine Species of Birds
}
\description{
The data consist of the egg boundary data of nine species of birds.
}
\usage{
data(eggs)
}
\details{
In the data set, there are four columns of variables: \code{Code}, \code{LatinName}, \code{x}, and \code{y}. 
\code{Code} saves the codes of individual eggs; 
\code{LatinName} saves the Latin names of the nine species of birds; 
\code{x} saves the \eqn{x} coordinates of the egg boundary in the Cartesian coordinate system (cm); 
and \code{y} saves the \eqn{y} coordinates of the egg boundary in the Cartesian coordinate system (cm). 
In \code{Code}, codes 1-9 represent \emph{Strix uralensis}, \emph{Dromaius novaehollandiae}, 
\emph{Turdus philomelos}, \emph{Gallus gallus}, \emph{Pandion haliaetus}, \emph{Uria aalge}, 
\emph{Uria lomvia}, \emph{Gallinago media}, and \emph{Aptenodytes patagonicus}, respectively. 
}
\references{
Narushin, V.G., Romanov, M.N., Griffin, D.K. (2021) Egg and math: introducing a universal formula for egg shape. 
\emph{Annals of the New York Academy of Sciences} 1505, 169\eqn{-}177. \doi{10.1111/nyas.14680}

Shi, P., Gielis, J., Niklas, K.J. (2022) Comparison of a universal (but complex) model for avian egg 
shape with a simpler model. \emph{Annals of the New York Academy of Sciences} 1514, 34\eqn{-}42. \doi{10.1111/nyas.14799}

Tian, F., Wang, Y., Sandhu, H.S., Gielis, J., Shi, P. (2020) Comparison of seed morphology of two ginkgo cultivars. 
\emph{Journal of Forestry Research} 31, 751\eqn{-}758. \doi{10.1007/s11676-018-0770-y}
}
\examples{
data(eggs)

uni.C <- sort( unique(eggs$Code) )
ind   <- 8
Data  <- eggs[eggs$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1 <- adjdata(x0, y0, ub.np=1000, times=1.2, len.pro=1/20)
x1   <- Res1$x
y1   <- Res1$y

dev.new()
plot( x1, y1, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
      xlab=expression(italic("x")), ylab=expression(italic("y")),
      pch=1, col=1 ) 

Res2 <- adjdata(x0, y0, ub.np=60, times=1, len.pro=1/2, index.sp=20)
x2   <- Res2$x
y2   <- Res2$y

Res3 <- adjdata(x0, y0, ub.np=60, times=1, len.pro=1/2, index.sp=100)
x3   <- Res3$x
y3   <- Res3$y

dev.new()
plot( x2, y2, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")),
      pch=1, col=4 ) 
points( x3, y3, col=2)

graphics.off()
}
