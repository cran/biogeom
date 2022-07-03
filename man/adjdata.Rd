\name{adjdata}
\encoding{UTF-8}
\alias{adjdata}
\title{
Boundary Data Adjustment of A Polygon
}
\description{
\code{adjdata} adjusts the data points in counterclockwise order based on the shortest distance method.
}
\usage{
adjdata(x, y, ub.np = 2000, times = 1.2, len.pro = 1/20, index.sp = 1)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of points on a polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of points on a polygon's boundary.}
  \item{ub.np}{the upper bound of the number of points eventually retained on the polygon's boundary.}
  \item{times}{the number of times \code{un.np} is initially retained as the number of points.}
  \item{len.pro}{the proportion of the distance between any two points to 
the maximum distance between the points on the polygon's boundary, which is used to 
determine whether the second point needs to be deleted.}
  \item{index.sp}{the index of the starting point of a group of indices that regularly 
divide the number of points on the polygon's boundary into \code{ub.np} parts.}
}
\details{
When \code{ub.np} > \code{length(x)}, \code{length(x)} points 
on the polygon's boundary are retained. 
The \code{\link{quantile}} function in package \pkg{stats} 
is used to carry out the regular division of data points. From the starting point, the second point is the 
one that has the shortest distance from the former. When the distance between the two points 
is larger than \code{len.pro} multiplied by the maximum distance between points on the polygon's boundary,   
the second point is deleted from the coordinates. 
Then, the third point that has the shortest distance from the first point is defined as the second point. 
If the distance between the first point and the second point is no more than \code{len.pro} multiplied 
by the maximum distance, the first and second points are recorded in a new matrix for the coordinates of the polygon, 
and the second point is defined as the first point in the old matrix for the coordinates of the polygon. 
The shortest distance method is then used to look for a third point that meets the requirement. 
}
\value{
  \item{x}{the \eqn{x} coordinates of points eventually retained on the polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of points eventually retained on the polygon's boundary.}
}
\references{
Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14862}

Shi, P., Ratkowsky, D.A., Li, Y., Zhang, L., Lin, S., Gielis, J. (2018) General leaf-area geometric 
formula exists for plants - Evidence from the simplified Gielis equation. 
\emph{Forests} 9, 714. \doi{10.3390/f9110714}

Su, J., Niklas, K.J., Huang, W., Yu, X., Yang, Y., Shi, P. (2019) 
Lamina shape does not correlate with lamina surface area: An analysis based on the simplified Gielis equation. 
\emph{Global Ecology and Conservation} 19, e00666. \doi{10.1016/j.gecco.2019.e00666}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\note{
The initial boundary data of a polygon can be obtained by running the M-file based on Matlab (version >= 2009a) 
developed by Shi et al. (2018) and Su et al. (2019) for a .bmp black and white image of the polygon. See references below.
}
\examples{
data(eggs)
uni.C1 <- sort( unique(eggs$Code) )
ind1   <- 2
Data1  <- eggs[eggs$Code==uni.C1[ind1], ]
x0     <- Data1$x
y0     <- Data1$y

Res1   <- adjdata(x0, y0, ub.np=2000, times=1.2, len.pro=1/20)
x1     <- Res1$x
y1     <- Res1$y

dev.new()
plot( x1, y1, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")),
      pch=1, col=1 ) 

Res2   <- adjdata(x0, y0, ub.np=40, times=1, len.pro=1/2, index.sp=20)
x2     <- Res2$x
y2     <- Res2$y

Res3   <- adjdata(x0, y0, ub.np=100, times=1, len.pro=1/2, index.sp=100)
x3     <- Res3$x
y3     <- Res3$y

dev.new()
plot( x2, y2, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")),
      pch=1, col=4 ) 
points( x3, y3, col=2)


data(starfish)

uni.C2 <- sort( unique(starfish$Code) )
ind2   <- 2
Data2  <- starfish[starfish$Code==uni.C2[ind2], ]
x4     <- Data2$x
y4     <- Data2$y

dev.new()
plot( x4, y4, asp=1, type="l", cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

Res4 <- adjdata(x4, y4, ub.np=500, times=1.2, len.pro=1/20)
x5   <- Res4$x
y5   <- Res4$y

dev.new()
plot( x5, y5, asp=1, type="l", cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

}