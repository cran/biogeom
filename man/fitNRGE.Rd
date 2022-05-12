\name{fitNRGE}
\alias{fitNRGE}
\title{
Parameter Estimation for the Narushin-Romanov-Griffin Equation
}
\description{
\code{fitNRGE} is used to estimate the parameters of the Narushin-Romanov-Griffin equation.
}
\usage{
fitNRGE(x, y, angle = NULL, x0 = NULL, y0 = NULL, ini.C = c(-1, 0.1, 0.5, 1), 
        strip.num = 2000, control = list(), fig.opt = TRUE, xlim = NULL,  
        ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of the edge of an egg's boundary.}
  \item{y}{the \eqn{y} coordinates of the edge of an egg's boundary.}
  \item{angle}{the angle between the major axis (i.e., the egg length axis) of the egg shape and the \eqn{x}-axis, which is estimated by SGE.}
  \item{x0}{the polar point's abscissa corresponding to \code{!is.null(angle)}, when using the SGE major axis approximation method.}
  \item{y0}{the polar point's ordinate corresponding to \code{!is.null(angle)}, when using the SGE major axis approximation method.}
  \item{ini.C}{the initial value(s) of parameter \eqn{C} in the Narushin-Romanov-Griffin equation.}  
  \item{strip.num}{the number of equidistant strips intersecting with 
the egg's boundary that are horizontally placed. See Shi et al. (2018, 2020) for details.}  
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{fig.opt}{an optional argument to draw the observed and predicted egg's boundaries.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Narushin-Romanov-Griffin curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Narushin-Romanov-Griffin curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Narushin-Romanov-Griffin curve.}
  \item{main}{the main title of the figure.}
}
\details{
The NRGE (see \code{\link{NRGE}}) has a complex model structure with four parameters (i.e., \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D}). 
Because three out of four parameters of NRGE have clear biological and geometric meanings (i.e., \eqn{A}, \eqn{B}, and \eqn{D}), 
their values could be estimated by means of numerical calculation. After obtaining the numerical values of the three parameters, 
the Nelder-Mead algorithm (Nelder and Mead, 1965) was used to estimate \eqn{C}. 
Because of the failure of the optimization method to estimate the major axis (i.e., the egg length axis) and model parameters of NRGE, 
it was difficult to define the egg length axis, although it is essential for calculating \eqn{A}, \eqn{B}, and \eqn{D}. 
For this reason, two methods were used to obtain the major axis: the maximum distance method, and 
the SGE major axis approximation method, where SGE represents a simplified Gielis equation (\code{\link{GE}} with \code{simpver = 1} and \code{m = 1}). 
In the first method, the straight line through two points forming the maximum distance 
on the egg's boundary is defined as the major axis. In the second method, the major axis predicted by SGE 
was directly used as the major axis of NRGE, because SGE balances the goodness-of-fit of the model and the bilateral symmetry of the curve. 
Because the direction from the egg base to the egg tip predicted by SGE is 
the reverse of that predicted by NRGE, the angle between the major axis of NRGE and the \eqn{x}-axis is equal to 
the sum of the estimated angle of the major axis using SGE and \eqn{\pi}. 
When \code{angle = NULL}, the maximum distance method is used; when \code{angle} is a numerical value, 
the SGE major axis approximation method is used, where \code{x0}, \code{y0}, and \code{angle} 
should be equal to the first three estimated parameters using the \code{\link{fitGE}} 
function with arguments \code{m = 1} and \code{simpver = 1}. Here, the numerical value of \code{angle} is not 
the angle between the major axis of NRGE and the \eqn{x}-axis, and instead it is the angle between the major axis of 
SGE and the \eqn{x}-axis. 
Once the major axis is established, the distance of the major axis can be calculated as the estimate of \eqn{A}. 
Using the maximum distance method, \eqn{A} equals the maximum distance. Using the SGE major axis approximation method, 
\eqn{A} may be slightly smaller than the true distance. After rotating the major axis to make it overlap with the \eqn{x}-axis, 
a large number of equidistant strips can be used (e.g., 2000) from the egg base to egg tip to intersect the egg's boundary. 
This methodology makes it easy to obtain the maximum egg width (i.e., \eqn{B}) and \eqn{D}. The residual sum of squares (RSS) 
between the observed and predicted \eqn{y} values can be minimized using an optimization method (Nelder and Mead, 1965) to estimate \eqn{C}. 
Despite the complex structure of NRGE (see \code{\link{NRGE}}), the optimization method for estimating the remaining parameter \eqn{C} 
becomes feasible after the other three parameters have been numerically estimated. Please see Shi et al. (2022) for details.
}
\value{
  \item{theta}{the angle between the major axis (i.e., the egg length axis) of the egg shape and the \eqn{x}-axis.}
  \item{x.obs}{the observed \eqn{x} coordinates.}
  \item{y.obs}{the observed \eqn{y} coordinates.}
  \item{y.pred}{the predicted \eqn{y} coordinates corresponding to the the observed \eqn{x} coordinates.}
  \item{par}{the estimates of the four model parameters in the Narushin-Romanov-Griffin equation.}
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
and it is equal to the sum of \eqn{\pi} and the given angle when using the SGE major axis approximation method 
(i.e., the angle between the SGE major axis and the \eqn{x}-axis). 
There are four estimated parameters in total for \code{par}.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{

Narushin, V.G., Romanov, M.N., Griffin, D.K. (2021) Egg and math: introducing a universal formula for egg shape. 
\emph{Annals of the New York Academy of Sciences} 1505, 169\eqn{-}177. \doi{10.1111/nyas.14680}

Nelder, J.A., Mead, R. (1965). A simplex algorithm for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Gielis, J., Niklas, K.J. (2022) Comparison of a universal (but complex) model for avian egg 
shape with a simpler model. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14799}

Shi, P., Niinemets, U., Hui, C., Niklas, K.J., Yu, X., Holscher, D. (2020) 
Leaf bilateral symmetry and the scaling of the perimeter vs. the surface area in 15 vine species. 
\emph{Forests} 11, 246. \doi{10.3390/f11020246}

Shi, P., Zheng, X., Ratkowsky, D.A., Li, Y., Wang, P., Cheng, L. (2018)  
A simple method for measuring the bilateral symmetry of leaves. 
\emph{Symmetry} 10, 118. \doi{10.3390/sym10040118}

Tian, F., Wang, Y., Sandhu, H.S., Gielis, J., Shi, P. (2020) Comparison of seed morphology of two ginkgo cultivars. 
\emph{Journal of Forestry Research} 31, 751\eqn{-}758. \doi{10.1007/s11676-018-0770-y}
}
\seealso{
\code{\link{curveNRGE}}, \code{\link{fitNRGE}}
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
  Res2 <- fitNRGE(x1, y1, angle=NULL, ini.C=c(-1, -0.1, seq(0.1, 1, by=0.05)), 
                  strip.num=2000, control=list(), fig.opt=TRUE) 
         
  dev.new()
  plot(Res2$x.obs, Res2$y.obs, asp=1, cex.lab=1.5, cex.axis=1.5, 
       xlab=expression(italic("x")), ylab=expression(italic("y")),
       type="l", col=4) 
  lines( Res2$x.obs, Res2$y.pre, col=2)
}

}

