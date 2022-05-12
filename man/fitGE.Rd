\name{fitGE}
\alias{fitGE}
\title{
Data-Fitting Function for the Gielis Equation
}
\description{
\code{fitGE} is used to estimate the parameters of the original (or twin) Gielis equation 
or one of its simplified versions.
}
\usage{
fitGE(expr, x, y, ini.val, m = 1, simpver = NULL, 
      nval = nval, control = list(), par.list = FALSE, 
      stand.fig = TRUE, angle = NULL, fig.opt = FALSE, np = 2000,
      xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{expr}{the original (or twin) Gielis equation or one of its simplified versions.}
  \item{x}{the \eqn{x} coordinates of a polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of a polygon's boundary.}
  \item{ini.val}{the list of initial values for the model parameters.}
  \item{m}{the given \eqn{m} value that determines the number of angles of the Gielis curve within \eqn{[0, 2\pi)}.}
  \item{simpver}{an optional argument to use the simplified version of the original (or twin) Gielis equation.}
  \item{nval}{the specified value for \eqn{n_{1}} or \eqn{n_{2}} or \eqn{n_{3}} in the simplified versions.}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{par.list}{the option of showing the list of parameters on the screen.}
  \item{stand.fig}{the option of drawing the observed and predicted polygons at the standard state 
(i.e., the polar point is located at (0, 0), and the major axis overlaps with the \eqn{x}-axis).}
  \item{angle}{the angle between the major axis and the \eqn{x}-axis, which can be defined by the user.}
  \item{fig.opt}{an optional argument of drawing the observed and predicted polygons at arbitrary angle 
between the major axis and the \eqn{x}-axis.}
  \item{np}{the number of data points on the predicted Gielis curve.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Gielis curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Gielis curve.}
  \item{unit}{the unit of the \eqn{x}-axis and the \eqn{y}-axis when showing the Gielis curve.}
  \item{main}{the main title of the figure.}
}
\details{
The arguments of \code{m}, \code{simpver}, and \code{nval} should correspond 
to \code{expr} (i.e., \code{GE} or \code{TGE}). Please note the differences in the simplified version number and the 
number of parameters between \code{GE} and \code{TGE}. 
The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization of minimizing 
the residual sum of squares (RSS) between the observed and predicted radii. 
The \code{\link{optim}} function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
When \code{angle = NULL}, the observed polygon will be shown at its initial angle in the scanned image; 
when \code{angle} is a numerical value (e.g., \eqn{\pi/4}) defined by the user, it indicates that the major axis 
is rotated by the amount (\eqn{\pi/4}) counterclockwise from the \eqn{x}-axis.
}
\value{
  \item{par}{the estimates of the model parameters.}
  \item{scan.length}{the observed length of the polygon.}
  \item{scan.width}{the observed width of the polygon.}
  \item{scan.area}{the observed area of the polygon.}
  \item{r.sq}{the coefficient of determination between the observed and predicted polar radii.}
  \item{RSS}{the residual sum of squares between the observed and predicted polar radii.}
  \item{sample.size}{the number of data points used in the data fitting.}
  \item{phi.stand.obs}{the polar angles at the standard state.}
  \item{phi.trans}{the transferred polar angles rotated as defined by the user.}
  \item{r.stand.obs}{the observed polar radii at the standard state.}
  \item{r.stand.pred}{the predicted polar radii at the standard state.}
  \item{x.stand.obs}{the observed \eqn{x} coordinates at the standard state.}
  \item{x.stand.pred}{the predicted \eqn{x} coordinates at the standard state.}
  \item{y.stand.obs}{the observed \eqn{y} coordinates at the standard state.}
  \item{y.stand.pred}{the predicted \eqn{y} coordinates at the standard state.} 
  \item{r.obs}{the observed polar radii at the transferred polar angles as defined by the user.}
  \item{r.pred}{the predicted polar radii at the transferred polar angles as defined by the user.}
  \item{x.obs}{the observed \eqn{x} coordinates at the transferred polar angles as defined by the user.}
  \item{x.pred}{the predicted \eqn{x} coordinates at the transferred polar angles as defined by the user.}
  \item{y.obs}{the observed \eqn{y} coordinates at the transferred polar angles as defined by the user.}
  \item{y.pred}{the predicted \eqn{y} coordinates at the transferred polar angles as defined by the user.}       
}
\note{
\code{simpver} in \code{GE} is different from that in \code{TGE}.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Gielis, J. (2003) A generic geometric transformation that unifies a wide range of natural 
and abstract shapes. \emph{American Journal of Botany} 90, 333\eqn{-}338. \doi{10.3732/ajb.90.3.333}

Li, Y., Quinn, B.K., Gielis, J., Li, Y., Shi, P. (2022) Evidence that supertriangles exist 
in nature from the vertical projections of \emph{Koelreuteria paniculata} fruit. \emph{Symmetry} 
14, 23. \doi{10.3390/sym14010023}

Nelder, J.A., Mead, R. (1965). A simplex algorithm for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Ratkowsky, D.A., Gielis, J. (2020) The generalized Gielis geometric
equation and its application. \emph{Symmetry} 12, 645. \doi{10.3390/sym12040645}

Shi, P., Xu, Q., Sandhu, H.S., Gielis, J., Ding, Y., Li, H., Dong, X. (2015) Comparison
of dwarf bamboos (\emph{Indocalamus} sp.) leaf parameters to determine relationship
between spatial density of plants and total leaf area per plant. \emph{Ecology and Evolution} 5,
4578-4589. \doi{10.1002/ece3.1728}
}
\seealso{
\code{\link{areaGE}}, \code{\link{curveGE}}, \code{\link{GE}}, \code{\link{TGE}}
}
\examples{
data(eggs)

uni.C <- sort( unique(eggs$Code) )
ind   <- 1
Data  <- eggs[eggs$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1 <- adjdata(x0, y0, ub.np=200, times=1.2, len.pro=1/20)
x1   <- Res1$x
y1   <- Res1$y
Res2 <- adjdata(x0, y0, ub.np=40, times=1, len.pro=1/2, index.sp=20)
x2   <- Res2$x
y2   <- Res2$y
Res3 <- adjdata(x0, y0, ub.np=100, times=1, len.pro=1/2, index.sp=100)
x3   <- Res3$x
y3   <- Res3$y

dev.new()
plot( x2, y2, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")),
      pch=1, col=4 ) 
points( x3, y3, col=2)

\donttest{
  x0.ini    <- mean( x1 )
  y0.ini    <- mean( y1 )
  theta.ini <- pi
  a.ini     <- sqrt(2) * max( y0.ini-min(y1), x0.ini-min(x1) ) 
  n1.ini <- c(5, 25)
  n2.ini <- c(15, 25)
  if(ind == 2){
    n1.ini <- c(0.5, 1)
    n2.ini <- c(6, 12)
  }
  ini.val <- list(x0.ini, y0.ini, theta.ini, a.ini, n1.ini, n2.ini)

  Res4 <- fitGE( GE, x=x1, y=y1, ini.val=ini.val, 
                 m=1, simpver=1, nval=1, unit="cm",  
                 par.list=FALSE, fig.opt=TRUE, angle=NULL, 
                 control=list(reltol=1e-20, maxit=20000), 
                 np=2000 )
  Res4$par
  sqrt(sum((Res4$y.stand.obs-Res4$y.stand.pred)^2)/Res4$sample.size)

  x4 <- Res4$x.stand.pred
  y4 <- Res4$y.stand.pred

  library(spatstat.geom)
  poly4 <- as.polygonal(owin(poly=list(x=x4, y=y4)))
  area(poly4)

  areaGE(GE, P = Res4$par[4:6], 
         m=1, simpver=1)
}

}