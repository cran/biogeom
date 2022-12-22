\name{fitovate}
\encoding{UTF-8}
\alias{fitovate}
\title{
Data-Fitting Function for the Ovate Leaf-Shape Equation
}
\description{
\code{fitovate} is used to estimate the parameters of a simplified performance equation.
}
\usage{
fitovate(expr, x, y, ini.val, 
         par.list = FALSE, stand.fig = TRUE, control = list(), 
         angle = NULL, fig.opt = FALSE, index.xmax = 3, np = 2000, 
         xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{expr}{the simplified version 1 of the performance equations.}
  \item{x}{the \eqn{x} coordinates of a polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of a polygon's boundary.}
  \item{ini.val}{the initial values of the simplified version 1 of a performance equation.}  
  \item{par.list}{an optional argument to show the list of parameters on the screen.}
  \item{stand.fig}{an optional argument to draw the observed and predicted polygons' boundaries at the standard state 
(i.e., the origin is located at (0, 0), and the major axis overlaps with the \eqn{x}-axis).}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{angle}{the angle between the major axis of the polygon and the \eqn{x}-axis, which can be defined by the user.}
  \item{fig.opt}{an optional argument to draw the observed and predicted polygons at an arbitrary angle 
between the major axis and the \eqn{x}-axis.}
  \item{index.xmax}{the specified index in parameters representing \eqn{x_{\mathrm{max}}}.}
  \item{np}{the number of data points on the predicted ovate leaf-shape curve.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the ovate leaf-shape curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the ovate leaf-shape curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the ovate leaf-shape curve.}
  \item{main}{the main title of the figure.}
}
\details{
\code{ini.val} is a list for seven parameters: three location parameters, and four model parameters, i.e., \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, 
\eqn{x_{\mathrm{max}}}, and \eqn{\delta}. This means that \code{expr} is limited to being the simplified version 1 (where \eqn{x_{\mathrm{min}} = 0})  
in \code{\link{MbetaE}}, \code{\link{MBriereE}}, and \code{\link{MLRFE}}. 
The initial values for the first three parameters in \code{ini.val} are location parameters, 
among which the first two are the planar coordinates of the transferred origin, 
and the third is the angle between the major axis of the polygon and the \eqn{x}-axis. 
The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization of minimizing 
the residual sum of squares (RSS) between the observed and predicted radii. 
The \code{\link{optim}} function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
When \code{angle = NULL}, the observed polygon will be shown at its initial angle in the scanned image; 
when \code{angle} is a numerical value (e.g., \eqn{\pi/4}) defined by the user, it indicates that the major axis 
is rotated by the amount (\eqn{\pi/4}) counterclockwise from the \eqn{x}-axis.
}
\value{
  \item{par}{the estimates of the model parameters.}
  \item{r.sq}{the coefficient of determination between the observed and predicted \eqn{y} values.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values.}
  \item{sample.size}{the number of data points on the polygon's boundary in the data fitting.}
  \item{scan.length}{the observed length of the polygon's boundary.}
  \item{scan.width}{the observed width of the polygon's boundary.}
  \item{scan.perimeter}{the observed perimeter of the polygon's boundary.}
  \item{scan.area}{the observed area of the polygon's boundary.}
  \item{pred.length}{the predicted length of the polygon's boundary.}
  \item{pred.width}{the predicted width of the polygon's boundary.}
  \item{pred.perimeter}{the predicted perimeter of the polygon's boundary.}
  \item{pred.area}{the predicted area of the polygon's boundary.}
  \item{x.stand.obs}{the observed \eqn{x} coordinates at the standard state.}
  \item{x.stand.pred}{the predicted \eqn{x} coordinates at the standard state.}
  \item{y.stand.obs}{the observed \eqn{y} coordinates at the standard state.}
  \item{y.stand.pred}{the predicted \eqn{y} coordinates at the standard state.} 
  \item{x.obs}{the observed \eqn{x} coordinates at the transferred angles defined by the user.}
  \item{x.pred}{the predicted \eqn{x} coordinates at the transferred angles defined by the user.}
  \item{y.obs}{the observed \eqn{y} coordinates at the transferred angles defined by the user.}
  \item{y.pred}{the predicted \eqn{y} coordinates at the transferred angles defined by the user.}       
}
\note{
There are seven parameters in total for the value of \code{par}. The transferred angle denotes the angle 
between the major axis and the \eqn{x}-axis.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Jin, J., Quinn, B.K., Shi, P. (2022) The modified Brière equation and its 
applications. \emph{Plants} 11, 1769. \doi{10.3390/plants11131769}

Li, Y., Zheng, Y., Ratkowsky, D.A., Wei, H., Shi, P. (2022) Application of an ovate leaf shape model 
to evaluate leaf bilateral asymmetry and calculate lamina centroid location. 
\emph{Frontiers in Plant Science} 12, 822907. \doi{10.3389/fpls.2021.822907}

Nelder, J.A., Mead, R. (1965) A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Yu, K., Niklas, K.J., Schrader, J., Song, Y., Zhu, R., Li, Y., Wei, H., Ratkowsky, D.A. (2021)  
A general model for describing the ovate leaf shape. \emph{Symmetry} 13, 1524. \doi{10.3390/sym13081524}
}
\seealso{
\code{\link{areaovate}}, \code{\link{curveovate}}, \code{\link{MbetaE}}, \code{\link{MBriereE}}, \code{\link{MLRFE}}
}
\examples{
data(Neocinnamomum)

uni.C <- sort( unique(Neocinnamomum$Code) )
ind   <- 2
Data  <- Neocinnamomum[Neocinnamomum$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1 <- adjdata(x0, y0, ub.np=200, len.pro=1/20)
x1   <- Res1$x
y1   <- Res1$y

dev.new()
plot( Res1$x, Res1$y, asp=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic("x")), ylab=expression(italic("y")) ) 

\donttest{
  x0.ini    <- min( x1 )
  y0.ini    <- min( y1 )
  theta.ini <- pi/4
  len.max   <- max( max(y1)-min(y1), max(x1)-min(x1) ) *2/sqrt(2)
  a.ini     <- c(0.1, 0.01, 0.001, 0.0001)
  m.ini     <- c(0.1, 0.5, 1, 2)
  x2.ini    <- len.max
  delta.ini <- c(0.5, 1)
  ini.val   <- list(x0.ini, y0.ini, theta.ini, a.ini, m.ini, x2.ini, delta.ini)

  Res2 <- fitovate(MBriereE, x=x1, y=y1, ini.val=ini.val,            
            par.list=FALSE, fig.opt=TRUE, angle=pi/6, 
            control=list(reltol=1e-20, maxit=20000), 
            np=2000, unit=NULL)
  Res2$RSS

  x0.ini    <- min( x1 )
  y0.ini    <- min( y1 )
  theta.ini <- pi/4
  len.max   <- max( max(y1)-min(y1), max(x1)-min(x1) ) *2/sqrt(2)
  yc.ini    <- len.max/3 
  xc.ini    <- 1/4*len.max
  x2.ini    <- len.max
  delta.ini <- c(0.5, seq(1, 5, by=5))
  ini.val   <- list(x0.ini, y0.ini, theta.ini, yc.ini, xc.ini, x2.ini, delta.ini)

  Res3 <- fitovate( MbetaE, x=x1, y=y1, ini.val=ini.val,            
                    par.list=TRUE, fig.opt=TRUE, angle=pi/3, 
                    control=list(reltol=1e-20, maxit=20000), 
                    np=2000, unit=NULL )
  Res3$RSS

  Res4 <- fitovate( MLRFE, x=x1, y=y1, ini.val=ini.val, 
                    unit=NULL, par.list=FALSE, fig.opt=TRUE, 
                    angle=NULL, control=list(reltol=1e-20, 
                    maxit=20000), np=2000)
  Res4$RSS
}

graphics.off()
}