\name{fitPE}
\encoding{UTF-8}
\alias{fitPE}
\title{
Data-Fitting Function for the Preston Equation
}
\description{
\code{fitPE} is used to estimate the parameters of the original Preston equation 
or one of its simplified versions.
}
\usage{
fitPE(x, y, ini.val, simpver = NULL, 
      control = list(), par.list = FALSE, 
      stand.fig = TRUE, angle = NULL, fig.opt = FALSE, np = 2000,
      xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of a polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of a polygon's boundary.}
  \item{ini.val}{the list of initial values for the model parameters.}
  \item{simpver}{an optional argument to use the simplified version of the original Preston equation.}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{par.list}{the option of showing the list of parameters on the screen.}
  \item{stand.fig}{the option of drawing the observed and predicted polygons at the standard state 
(i.e., the polar point is located at (0, 0), and the major axis overlaps with the \eqn{x}-axis).}
  \item{angle}{the angle between the major axis and the \eqn{x}-axis, which can be defined by the user.}
  \item{fig.opt}{an optional argument of drawing the observed and predicted polygons at arbitrary angle 
between the major axis and the \eqn{x}-axis.}
  \item{np}{the number of data points on the predicted Preston curve.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Preston curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Preston curve.}
  \item{unit}{the unit of the \eqn{x}-axis and the \eqn{y}-axis when showing the Preston curve.}
  \item{main}{the main title of the figure.}
}
\details{
The \code{simpver} argument should correspond to \code{PE}. Here, the major axis is a straight line through 
the midpoint of the egg length, and is perpendicular to the egg length axis.  
The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization of minimizing 
the residual sum of squares (RSS) between the observed and predicted distances from the origin. 
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
  \item{r.sq}{the coefficient of determination between the observed and predicted distances 
of the points on the Preston curve from the origin.}
  \item{RSS}{the residual sum of squares between the observed and predicted distances 
of the points on the Preston curve from the origin.}
  \item{sample.size}{the number of data points used in the data fitting.}
  \item{zeta.stand.obs}{the \code{zeta} angles at the standard state.}
  \item{zeta.trans}{the transferred angles rotated as defined by the user.}
  \item{r.stand.obs}{the observed distances of the points on the Preston curve from the origin at the standard state.}
  \item{r.stand.pred}{the predicted distances of the points on the Preston curve from the origin at the standard state.}
  \item{x.stand.obs}{the observed \eqn{x} coordinates of the points on the Preston curve at the standard state.}
  \item{x.stand.pred}{the predicted \eqn{x} coordinates of the points on the Preston curve at the standard state.}
  \item{y.stand.obs}{the observed \eqn{y} coordinates of the points on the Preston curve at the standard state.}
  \item{y.stand.pred}{the predicted \eqn{y} coordinates of the points on the Preston curve at the standard state.} 
  \item{r.obs}{the observed distances of the points on the Preston curve from the origin 
at the transferred polar angles as defined by the user.}
  \item{r.pred}{the predicted distances of the points on the Preston curve from the origin 
at the transferred polar angles as defined by the user.}
  \item{x.obs}{the observed \eqn{x} coordinates of the points on the Preston curve 
at the transferred polar angles as defined by the user.}
  \item{x.pred}{the predicted \eqn{x} coordinates of the points on the Preston curve 
at the transferred polar angles as defined by the user.}
  \item{y.obs}{the observed \eqn{y} coordinates of the points on the Preston curve 
at the transferred polar angles as defined by the user.}
  \item{y.pred}{the predicted \eqn{y} coordinates of the points on the Preston curve 
at the transferred polar angles as defined by the user.}       
}
\note{
When \eqn{c_{1}} is a positive number, the egg tip generated by the Preston equation is down; 
when \eqn{c_{1}} is a negative number, the egg tip is up. To conveniently quantify the angle between 
the major axis and the \eqn{x}-axis, we define \eqn{c_{1}} to be a non-negative number. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Biggins, J.D., Thompson, J.E., Birkhead, T.R. (2018) Accurately quantifying 
the shape of birds' eggs. \emph{Ecology and Evolution} 8, 9728\eqn{-}9738. \doi{10.1002/ece3.4412}

Nelder, J.A., Mead, R. (1965). A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14862}
}
\seealso{
\code{\link{curvePE}}, \code{\link{PE}}
}
\examples{
data(eggs)

uni.C <- sort( unique(eggs$Code) )
ind   <- 9
Data  <- eggs[eggs$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

Res1 <- adjdata(x0, y0, ub.np=2000, times=1.2, len.pro=1/20)
x1   <- Res1$x
y1   <- Res1$y

dev.new()
plot( x1, y1, asp=1, cex.lab=1.5, cex.axis=1.5,  type="l", col=4,  
      xlab=expression(italic("x")), ylab=expression(italic("y")) )

\donttest{
  x0.ini    <- mean( x1 )
  y0.ini    <- mean( y1 )
  theta.ini <- -pi/4
  if(ind==1) theta.ini <- 0
  if(ind==7 | ind==9) theta.ini <- 4/3*pi    
  if(ind==8) theta.ini <- pi    
  a.ini     <- max(c(max(y1)-min(y1), max(x1)-min(x1)))/2
  b.ini     <- a.ini*3/4
  c1.ini    <- log(0.25)
  c2.ini    <- c(-0.1, -0.01, -0.001)  

  simpver   <- 1
  ini.val   <- list(x0.ini, y0.ini, theta.ini, a.ini, b.ini, c1.ini, c2.ini)

  res0 <- fitPE( x=x1, y=y1, ini.val=ini.val, 
                 simpver=simpver, unit="cm", par.list=FALSE, 
                 stand.fig=FALSE, angle=NULL, fig.opt=FALSE, 
                 control=list(reltol=1e-30, maxit=50000), np=2000 )

  n.loop <- 10
  Show   <- "FALSE"
  for(i in 1:n.loop){
    ini.val <- res0$par
    if(i==n.loop) Show <- "TRUE"
    print(paste(i, "/", n.loop, sep=""))
    res0 <- fitPE( x=x1, y=y1, ini.val=ini.val, 
                   simpver=simpver, unit="cm", par.list=FALSE, 
                   stand.fig=Show, angle=pi/2, fig.opt=Show,  
                   control=list(reltol=1e-30, maxit=50000), np=2000 )    
  }

  # The numerical values of the location and model parameters
  res0$par
  # The root-mean-square error (RMSE) between the observed and predicted distances from the origin
  sqrt(res0$RSS/res0$sample.size)
  # The root-mean-square error (RMSE) between the observed and predicted x coordinates
  sqrt(sum((res0$x.stand.obs-res0$x.stand.pred)^2)/length(res0$x.stand.obs))
}

}