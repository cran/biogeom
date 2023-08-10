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
lmPE(x, y, simpver = NULL, dev.angle = NULL, weights = NULL, fig.opt = TRUE, 
     prog.opt = TRUE, xlim = NULL, ylim = NULL, unit = NULL, main = NULL, 
     extr.method = "Shi")
}
\arguments{
  \item{x}{the \eqn{x} coordinates of the edge of an egg's boundary.}
  \item{y}{the \eqn{y} coordinates of the edge of an egg's boundary.}
  \item{simpver}{an optional argument to use the simplified version of the original Todd-Smart equation.}
  \item{dev.angle}{the angle of deviation for the axis associated with the maximum distance between two points 
on an egg's profile from the mid-line of the egg's profile.}
  \item{weights}{the weights for the multiple linear regression.}
  \item{fig.opt}{an optional argument to draw the observed and predicted egg's boundaries.}
  \item{prog.opt}{an optional argument to show the running progress for different values of \code{dev.angle}, when \code{dev.angle} is not \code{NULL}.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Todd-Smart curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Todd-Smart curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Todd-Smart curve.}
  \item{main}{the main title of the figure.}
  \item{extr.method}{an optional argument to fit the planar coordinate data of an egg's profile extracted using different methods.}
}
\details{
There are two methods to obtain the major axis (i.e., the mid-line) of an egg's profile: the maximum distance method, and 
the longest axis adjustment method. In the first method, the straight line through two points forming the maximum distance 
on the egg's boundary is defined as the major axis. In the second method, we assume that there is an angle 
of deviation for the longest axis (i.e., the axis associated with the maximum distance between two points on an egg's profile) 
from the mid-line of the egg's profile. That is to say, 
the mid-line of an egg's profile is not the axis associated with the maximum distance between two points on the egg's profile. 
When \code{dev.angle = NULL}, the maximum distance method is used; when \code{dev.angle} is a numerical value or a numerical vector, 
the longest axis adjustment method is used. Here, the numerical value of \code{dev.angle} is not 
the angle of deviation for the major axis of an egg's profile from the \eqn{x}-axis, and instead it is the angle of deviation  
for the longest axis (associated with the maximum distance between two points on the egg's profile) from the mid-line of the egg's profile. 
It is better to take the option of \code{extr.method = "Shi"} for correctly fitting the planar coordinate data of an egg's profile extracted 
using the protocols proposed by Shi et al. (2015, 2018) (and also see Su et al. (2019)), 
while it is better to take the option of \code{extr.method = "Biggins"} for correctly fitting the planar coordinate data 
of an egg's profile extracted using the protocols proposed by Biggins et al. (2018). For the planar coordinate data extracted using 
the protocols of Biggins et al. (2018), there are fewer data points on the two ends of the mid-line than other parts of the egg's profile, 
which means that the range of the observed \eqn{x} values might be smaller than the actual egg's length. A group of equidistant \eqn{x} 
values are set along the mid-line, and each \eqn{x} value corresponds to two \eqn{y} values that are respectively located at the upper 
and lower sides of the egg's profile. Because of the difference in the curvature for differnt parts of the egg's profile, 
the equidistant \eqn{x} values cannot render the extracted data points on the egg's profile to be regular. For the planar coordinate data 
extracted using the protocols of Shi et al. (2015, 2018), the data points are more regularly distributed on the egg's profile (perimeter) 
than those of Biggins et al. (2018), although the \eqn{x} values of the data points along the mid-line are not equidistant. 
}
\value{
  \item{lm.tse}{the fitted results of the multiple linear regression.}
  \item{par}{the estimates of the four model parameters in the Todd-Smart equation.}
  \item{theta}{the angle between the longest axis of an egg's profile (i.e., 
the axis associated with the maximum distance between two points on the egg's profile) and the \eqn{x}-axis.}
  \item{epsilon}{the optimal angle of deviation for the longest axis (associated with the maximum distance between two points on an egg's profile) 
from the mid-line of the egg's profile, when \code{dev.angle} is not \code{NULL}.}
  \item{RSS.vector}{the vector of residual sum of squares corresponding to \code{dev.angle}, when \code{dev.angle} is not \code{NULL}.}
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
  \item{RSS.scaled}{the residual sum of squares between the observed and predicted \eqn{y} values for a scaled egg's profile whose length equals 2.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values.}
  \item{sample.size}{the number of data points used in the numerical calculation.}
  \item{RMSE.scaled}{the root-mean-square errors between the observed and predicted \eqn{y} values for a scaled egg's profile whose length equals 2.}      
  \item{RMSE}{the root-mean-square errors between the observed and predicted \eqn{y} values.}       
}
\note{
\code{theta} is the calculated angle between the longest axis (i.e., the axis associated with the maximum distance 
between two points on an egg's profile) and the \eqn{x}-axis, and \code{epsilon} is the calculated angle of deviation for the longest 
axis from the mid-line of the egg's profile. This means that the angle between the mid-line and the \eqn{x}-axis is equal to 
\code{theta} + \code{epsilon}. 
Here, \code{RSS}, and \code{RMSE} are for the observed and predicted \eqn{y} coordinates of the egg's profile, not for those 
when the egg's length is scaled to 2. There are two figures when \code{fig.opt = TRUE}: (i) the observed 
and predicted egg's boundaries when the egg's length is scaled to 2, and (ii) the observed 
and predicted egg's boundaries at their actual scales. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Biggins, J.D., Montgomeries, R.M., Thompson, J.E., Birkhead, T.R. (2022) 
Preston’s universal formula for avian egg shape. \emph{Ornithology} 
139, ukac028. \doi{10.1093/ornithology/ukac028}

Biggins, J.D., Thompson, J.E., Birkhead, T.R. (2018) Accurately quantifying 
the shape of birds' eggs. \emph{Ecology and Evolution} 8, 9728\eqn{-}9738. \doi{10.1002/ece3.4412}

Nelder, J.A., Mead, R. (1965). A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Huang, J., Hui, C., Grissino-Mayer, H.D., Tardif, J., Zhai, L., Wang, F., Li, B. (2015) 
Capturing spiral radial growth of conifers using the superellipse to model tree-ring geometric shape. 
\emph{Frontiers in Plant Science} 6, 856. \doi{10.3389/fpls.2015.00856}

Shi, P., Niinemets, Ü., Hui, C., Niklas, K.J., Yu, X., Hölscher, D. (2020) 
Leaf bilateral symmetry and the scaling of the perimeter vs. the surface area in 15 vine species. 
\emph{Forests} 11, 246. \doi{10.3390/f11020246}

Shi, P., Ratkowsky, D.A., Li, Y., Zhang, L., Lin, S., Gielis, J. (2018) General leaf-area geometric 
formula exists for plants - Evidence from the simplified Gielis equation. 
\emph{Forests} 9, 714. \doi{10.3390/f9110714}

Shi, P., Wang, L., Quinn, B.K., Gielis, J. (2023) A new program to estimate the parameters 
of Preston's equation, a general formula for describing the egg shape of birds. \emph{Symmetry} 
15, 231. \doi{10.3390/sym15010231}

Su, J., Niklas, K.J., Huang, W., Yu, X., Yang, Y., Shi, P. (2019) 
Lamina shape does not correlate with lamina surface area: An analysis based on the simplified Gielis equation. 
\emph{Global Ecology and Conservation} 19, e00666. \doi{10.1016/j.gecco.2019.e00666}

Todd, P.H., Smart, I.H.M. (1984) The shape of birds' eggs. \emph{Journal of Theoretical Biology} 
106, 239\eqn{-}243. \doi{10.1016/0022-5193(84)90021-3}
}
\seealso{
\code{\link{curveEPE}}, \code{\link{fitEPE}}, \code{\link{PE}}, \code{\link{TSE}} 
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
  Res2 <- lmPE(x1, y1, simpver=NULL, dev.angle=NULL, unit="cm")
  summary( Res2$lm.tse )
  Res2$RMSE.scaled / 2

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

  Res3 <- lmPE(x1, y1, simpver=NULL, dev.angle=seq(-0.05, 0.05, by=0.0001), unit="cm")
  summary( Res3$lm.tse )
  Res3$epsilon
  Res3$RMSE.scaled / 2

  Res4 <- lmPE(x1, y1, simpver=1, dev.angle=NULL, unit="cm")
  summary( Res4$lm.tse )

}

graphics.off()
}

