\name{lmPE}
\encoding{UTF-8}
\alias{lmPE}
\title{
Parameter Estimation for the Todd-Smart Equation
}
\description{
\code{lmPE} is used to estimate the parameters of the Todd-Smart equation using the multiple linear regression.
}
\usage{
lmPE(x, y, simpver = NULL, angle = NULL, x0 = NULL, y0 = NULL, strip.num = 2000, 
     weights = NULL, fig.opt = TRUE, xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of the edge of an egg's boundary.}
  \item{y}{the \eqn{y} coordinates of the edge of an egg's boundary.}
  \item{simpver}{an optional argument to use the simplified version of the original Todd-Smart equation.}
  \item{angle}{the angle between the major axis (i.e., the axis perpendicular to the egg length axis) of the egg shape 
and the \eqn{x}-axis, which is estimated by the \code{\link{fitPE}} function.}
  \item{x0}{the polar point's abscissa corresponding to \code{!is.null(angle)}, when using the PE major axis approximation method.}
  \item{y0}{the polar point's ordinate corresponding to \code{!is.null(angle)}, when using the PE major axis approximation method.}
  \item{strip.num}{the number of equidistant strips intersecting with 
the egg's boundary that are horizontally placed. See Shi et al. (2018, 2020) for details.}  
  \item{weights}{the weights for the multiple linear regression.}
  \item{fig.opt}{an optional argument to draw the observed and predicted egg's boundaries.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Todd-Smart curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Todd-Smart curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Todd-Smart curve.}
  \item{main}{the main title of the figure.}
}
\details{
There are two methods to obtain the major axis (i.e., the egg length axis): the maximum distance method, and 
the PE major axis approximation method. In the first method, the straight line through two points forming the maximum distance 
on the egg's boundary is defined as the major axis. In the second method, the major axis predicted by using the \code{\link{fitPE}} function 
was directly used as the major axis of TSE, because the optimimzation method based on the \code{\link{fitPE}} function balances 
the goodness-of-fit of the model and the bilateral symmetry of the curve. 
Because the direction from the egg base to the egg tip predicted by the \code{\link{fitPE}} function is 
perpendicular to that predicted by TSE, the angle between the major axis of TSE and the \eqn{x}-axis is equal to 
the difference between the estimated angle of the major axis using the \code{\link{fitPE}} function and \eqn{\pi/2}. 
However, there is no need to subtract \eqn{\pi/2} in the \code{angule} argument. 
When \code{angle = NULL}, the maximum distance method is used; when \code{angle} is a numerical value, 
the PE major axis approximation method is used, where \code{x0}, \code{y0}, and \code{angle} 
should be equal to the first three estimated parameters using the \code{\link{fitPE}} 
function. Here, the numerical value of \code{angle} is not 
the angle between the major axis of TSE and the \eqn{x}-axis, and instead it is the angle between the major axis of 
PE and the \eqn{x}-axis. 
}
\value{
  \item{lm.tse}{the fitted results of the multiple linear regression.}
  \item{par}{the estimates of the four model parameters in the Todd-Smart equation.}
  \item{theta}{the angle between the major axis (i.e., the egg length axis) of the egg shape and the \eqn{x}-axis.}
  \item{x.obs}{the observed \eqn{x} coordinates.}
  \item{y.obs}{the observed \eqn{y} coordinates.}
  \item{y.pred}{the predicted \eqn{y} coordinates corresponding to the the observed \eqn{x} coordinates.}
  \item{x.stand.obs}{the observed \eqn{x} coordinates when the egg length is fixed to be 2 ranging from -1 to 1.}
  \item{y.stand.obs}{the observed \eqn{y} coordinates when the egg length is fixed to be 2 ranging from -1 to 1.}
  \item{y.stand.pred}{the predicted \eqn{y} coordinates corresponding to the the observed \eqn{x} coordinates, 
when the egg length is fixed to be 2 ranging from -1 to 1.}
  \item{scan.length}{the length of the egg's boundary. The default is the maximum distance between two points on 
the egg's boundary.}
  \item{scan.width}{the maximum width of the egg's boundary.}
  \item{scan.area}{the area of the egg's boundary.}
  \item{scan.perimeter}{the perimeter of the egg's boundary based on all data points on the egg's boundary.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values.}
  \item{sample.size}{the number of data points used in the numerical calculation.}
  \item{RMSE}{the root-mean-square errors between the observed and predicted \eqn{y} values.}       
}
\note{
\code{theta} is the calculated angle between the egg length axis and the \eqn{x}-axis when using the maximum distance method; 
and it is equal to the difference between the given (fitted) angle and \eqn{\pi/2} when using the PE major axis approximation method 
(i.e., the angle between the PE major axis and the \eqn{x}-axis). Here, \code{RSS}, and \code{RMSE} are for 
the observed and predicted \eqn{y} coordinates of the egg shape, not for those 
when the egg length is fixed to be 2 ranging from -1 to 1. There are two figures when \code{fig.opt = TRUE}: (i) the observed 
and predicted egg boundaries when the egg length is fixed to be 2 ranging from -1 to 1, and (ii) the observed 
and predicted egg boundaries at their actual scales. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Biggins, J.D., Montgomeries, R.M., Thompson, J.E., Birkhead, T.R. (2022) 
Preston’s universal formula for avian egg shape. \emph{Ornithology} 
In press. \doi{10.1093/ornithology/ukac028}

Biggins, J.D., Thompson, J.E., Birkhead, T.R. (2018) Accurately quantifying 
the shape of birds' eggs. \emph{Ecology and Evolution} 8, 9728\eqn{-}9738. \doi{10.1002/ece3.4412}

Nelder, J.A., Mead, R. (1965). A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14862}

Shi, P., Niinemets, Ü., Hui, C., Niklas, K.J., Yu, X., Hölscher, D. (2020) 
Leaf bilateral symmetry and the scaling of the perimeter vs. the surface area in 15 vine species. 
\emph{Forests} 11, 246. \doi{10.3390/f11020246}

Shi, P., Zheng, X., Ratkowsky, D.A., Li, Y., Wang, P., Cheng, L. (2018)  
A simple method for measuring the bilateral symmetry of leaves. 
\emph{Symmetry} 10, 118. \doi{10.3390/sym10040118}

Todd, P.H., Smart, I.H.M. (1984) The shape of birds' eggs. \emph{Journal of Theoretical Biology} 
106, 239\eqn{-}243. \doi{10.1016/0022-5193(84)90021-3}
}
\seealso{
\code{\link{curvePE}}, \code{\link{fitPE}}, \code{\link{PE}}, \code{\link{TSE}} 
}
\examples{
data(eggs)

uni.C <- sort( unique(eggs$Code) )
ind   <- 8
Data  <- eggs[eggs$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1  <- adjdata(x0, y0, ub.np=3000, len.pro=1/20)
x1    <- Res1$x
y1    <- Res1$y

dev.new()
plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) ) 

\donttest{
  Res2 <- lmPE(x1, y1, simpver=NULL, angle=NULL, unit="cm")
  summary( Res2$lm.tse )

  if(FALSE){
    dev.new()
    xg1 <- seq(-1, 1, len=1000)
    yg1 <- TSE(P=Res2$par, x=xg1, simpver=NULL)
    xg2 <- seq(1, -1, len=1000)
    yg2 <- -TSE(P=Res2$par, x=xg2, simpver=NULL)
    plot(xg1, yg1, asp=1, type="l", col=2, ylim=c(-1,1), cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic(x)), ylab=expression(italic(y)))
    lines(xg2, yg2, col=4)

    dev.new()
    plot(Res2$x.obs, Res2$y.obs, asp=1, cex.lab=1.5, cex.axis=1.5,  
      xlab=expression(italic(x)), ylab=expression(italic(y)), type="l")
    lines(Res2$x.obs, Res2$y.pred, col=2)

    dev.new()
    plot(Res2$x.stand.obs, Res2$y.stand.obs, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic(x)), ylab=expression(italic(y)), type="l")
    lines(Res2$x.stand.obs, Res2$y.stand.pred, col=2)
  }

  Res3 <- lmPE(x1, y1, simpver=1, angle=NULL, unit="cm")
  summary( Res3$lm.tse )

}

}

