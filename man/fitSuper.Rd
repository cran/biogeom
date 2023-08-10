\name{fitSuper}
\encoding{UTF-8}
\alias{fitSuper}
\title{
Data-Fitting Function for the Superellipse Equation
}
\description{
\code{fitSuper} is used to estimate the parameters of the superellipse equation.
}
\usage{
fitSuper(x, y, ini.val, control = list(), par.list = FALSE, 
         stand.fig = TRUE, angle = NULL, fig.opt = FALSE, np = 2000,
         xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of a polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of a polygon's boundary.}
  \item{ini.val}{the list of initial values for the model parameters.}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{par.list}{the option of showing the list of parameters on the screen.}
  \item{stand.fig}{the option of drawing the observed and predicted polygons at the standard state 
(i.e., the polar point is located at (0, 0), and the major axis overlaps with the \eqn{x}-axis).}
  \item{angle}{the angle between the major axis and the \eqn{x}-axis, which can be defined by the user.}
  \item{fig.opt}{an optional argument of drawing the observed and predicted polygons at arbitrary angle 
between the major axis and the \eqn{x}-axis.}
  \item{np}{the number of data points on the predicted superellipse curve.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the superellipse curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the superellipse curve.}
  \item{unit}{the unit of the \eqn{x}-axis and the \eqn{y}-axis when showing the superellipse curve.}
  \item{main}{the main title of the figure.}
}
\details{
The superellipse equation has its mathematical expression: 
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\varphi\right)\right|^{n}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\varphi\right)\right|^{n}\right)^{-\frac{1}{n}},}
where \eqn{r} represents the polar radius at the polar angle \eqn{\varphi}. 
\eqn{a} represents the semi-major axis of the superellipse; 
\eqn{k} is the ratio of the semi-minor axis to the semi-major axis of the superellipse; and \eqn{n} determines 
the curvature of the superellipse curve. This function is basically equal to the \code{\link{fitGE}} 
function with the arguments \code{m = 4} and \code{simpver = 9}. 
However, this function can make the estimated value of the parameter \eqn{k} be smaller than or equal to 1. 
Apart from the above three model parameters, there are three additional location parameters, i.e., the 
planar coordinates of the polar point (\eqn{x_{0}} and \eqn{y_{0}}), and the angle between the major axis 
of the superellipse and the \eqn{x}-axis (\eqn{\theta}). The input order of \code{ini.val} is \eqn{x_{0}}, \eqn{y_{0}}, 
\eqn{\theta}, \eqn{a}, \eqn{k}, and \eqn{n}. The fitted parameters will be shown after running this function 
in the same order. The Nelder-Mead algorithm (Nelder and Mead, 1965) 
is used to carry out the optimization of minimizing the residual sum of squares (RSS) between the observed and predicted radii. 
The \code{\link{optim}} function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
When \code{angle = NULL}, the observed polygon will be shown at its initial angle in the scanned image; 
when \code{angle} is a numerical value (e.g., \eqn{\pi/4}) defined by the user, it indicates that the major axis 
is rotated by the amount (\eqn{\pi/4}) counterclockwise from the \eqn{x}-axis.
}
\value{
  \item{par}{the estimated values of the parameters in the superellipse equation.}
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
The output values of running this function can be used as those of running the \code{\link{fitGE}} function. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Huang, J., Hui, C., Grissino-Mayer, H.D., Tardif, J., Zhai, L., Wang, F., Li, B. (2015) 
Capturing spiral radial growth of conifers using the superellipse to model tree-ring geometric shape. 
\emph{Frontiers in Plant Science} 6, 856. \doi{10.3389/fpls.2015.00856}

}
\seealso{
\code{\link{areaGE}}, \code{\link{curveGE}}, \code{\link{fitGE}}, \code{\link{GE}} 
}
\examples{
data(whitespruce)

uni.C <- sort( unique(whitespruce$Code) )  
Data  <- whitespruce[whitespruce$Code==uni.C[12], ]
x0    <- Data$x
y0    <- Data$y
Res1  <- adjdata(x0, y0, ub.np=200, len.pro=1/20)
x1    <- Res1$x
y1    <- Res1$y

plot( x1, y1, asp=1, cex.lab=1.5, cex.axis=1.5, type="l", 
      col="grey73", lwd=2,
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

\donttest{
  x0.ini    <- mean( x1 )
  y0.ini    <- mean( y1 )
  theta.ini <- c(0, pi/4, pi/2)
  a.ini     <- mean(c(max(x1)-min(x1), max(y1)-min(y1)))/2
  k.ini     <- 1
  n.ini     <- c(1.5, 2, 2.5)
  ini.val   <- list( x0.ini, y0.ini, theta.ini, a.ini, k.ini, n.ini )
  Res2      <- fitSuper(x=x1, y=y1, ini.val=ini.val, unit="cm", par.list=FALSE,
                        stand.fig=FALSE, angle=NULL, fig.opt=TRUE,  
                        control=list(reltol=1e-20, maxit=20000), np=2000)
  Res2$par
  Res2$r.sq
}

graphics.off()
}