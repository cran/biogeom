\name{fitETE}
\encoding{UTF-8}
\alias{fitETE}
\title{
Data-Fitting Function for the Explicit Troscianko Equation
}
\description{
\code{fitETE} is used to estimate the parameters of the explicit Troscianko equation.
}
\usage{
fitETE(x, y, ini.val, control = list(), par.list = FALSE, 
      stand.fig = TRUE, angle = NULL, fig.opt = FALSE, np = 2000,
      xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of an egg's profile.}
  \item{y}{the \eqn{y} coordinates of an egg's profile.}
  \item{ini.val}{the list of initial values for the model parameters.}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{par.list}{the option of showing the list of parameters on the screen.}
  \item{stand.fig}{the option of drawing the observed and predicted profiles of an egg at the standard state 
(i.e., the egg's centre is located at (0, 0), and the mid-line is aligned to the \eqn{x}-axis).}
  \item{angle}{the angle between the mid-line and the \eqn{x}-axis, which can be defined by the user.}
  \item{fig.opt}{an optional argument of drawing the observed and predicted profiles of an egg 
at arbitrary angle between the major axis and the \eqn{x}-axis.}
  \item{np}{the number of data points on the predicted explicit Troscianko curve.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Troscianko curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Troscianko curve.}
  \item{unit}{the unit of the \eqn{x}-axis and the \eqn{y}-axis when showing the Troscianko curve.}
  \item{main}{the main title of the figure.}
}
\details{
Here, the major axis (i.e., the mid-line of an egg's profile) is the straight line trhough the two ends of 
the egg's length. The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization 
of minimizing the residual sum of squares (RSS) between the observed and predicted \eqn{y} values. 
The \code{\link{optim}} function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
When \code{angle = NULL}, the observed egg's profile will be shown at its initial angle in the scanned image; 
when \code{angle} is a numerical value (e.g., \eqn{\pi/4}) defined by the user, it indicates that the major axis 
is rotated by the amount (\eqn{\pi/4}) counterclockwise from the \eqn{x}-axis.
}
\value{
  \item{par}{the estimates of the model parameters.}
  \item{scan.length}{the observed length of the egg's profile.}
  \item{scan.width}{the observed width of the egg's profile.}
  \item{scan.area}{the observed area of the egg's profile.}
  \item{scan.perimeter}{the observed perimeter of the egg's profile.}
  \item{r.sq}{the coefficient of determination between the observed and predicted \eqn{y} values 
on the Troscianko curve.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values 
on the Troscianko curve.}
  \item{sample.size}{the number of data points used in the data fitting.}
  \item{x.stand.obs}{the observed \eqn{x} coordinates of the points 
on the Troscianko curve at the standard state.}
  \item{y.stand.obs}{the observed \eqn{y} coordinates of the points 
on the Troscianko curve at the standard state.}
  \item{y.stand.pred}{the predicted \eqn{y} coordinates of the points 
on the Troscianko curve at the standard state.} 
  \item{x.obs}{the observed \eqn{x} coordinates of the points on the Troscianko curve 
at the transferred polar angles as defined by the user.}
  \item{y.obs}{the observed \eqn{y} coordinates of the points on the Troscianko curve 
at the transferred polar angles as defined by the user.}
  \item{y.pred}{the predicted \eqn{y} coordinates of the points on the Troscianko curve 
at the transferred polar angles as defined by the user.}       
}
\note{
In the outputs, there are no \code{x.stand.pred} and \code{x.pred}, because \code{y.stand.obs} and 
\code{y.stand.pred} share the same \eqn{x} values (i.e., \code{x.stand.obs}), and \code{y.obs} and 
\code{y.pred} share the same \eqn{x} values (i.e., \code{x.obs}).
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

Nelder, J.A., Mead, R. (1965) A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Wang, L., Quinn, B.K., Gielis, J. (2023) A new program to estimate the parameters 
of Preston's equation, a general formula for describing the egg shape of birds. \emph{Symmetry} 
15, 231. \doi{10.3390/sym15010231}

Troscianko, J. (2014). A simple tool for calculating egg shape, volume and surface area from digital images. 
\emph{Ibis}, 156, 874\eqn{-}878. \doi{10.1111/ibi.12177} 
}
\seealso{
\code{\link{curveETE}}, \code{\link{TE}}, \code{\link{lmTE}}
}
\examples{
data(eggs)

uni.C <- sort( unique(eggs$Code) )
ind   <- 8
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
  res1      <- lmTE( x1, y1, unit="cm", fig.opt=FALSE )
  
  if(FALSE){
    P0  <- c(res1$scan.length/2, res1$par)
    xx  <- seq(-res1$scan.length/2, res1$scan.length/2, len=2000)
    yy1 <- ETE(P0, xx)
    yy2 <- -ETE(P0, xx)
    dev.new()
    plot( xx, yy1, cex.lab=1.5, cex.axis=1.5, asp=1, col=2, 
          ylim=c(-res1$scan.length/2, res1$scan.length/2),
          type="l", xlab=expression(x), ylab=expression(y) )
    lines( xx, yy2, col=4 )
  }

  x0.ini    <- mean( x1 )
  y0.ini    <- mean( y1 )
  theta.ini <- res1$theta
  a.ini     <- res1$scan.length / 2
  alpha0.ini <- res1$par[1] 
  alpha1.ini <- res1$par[2]
  alpha2.ini <- res1$par[3]

  ini.val <- list(x0.ini, y0.ini, theta.ini, a.ini, alpha0.ini, alpha1.ini, alpha2.ini)

  res0 <- fitETE( x=x1, y=y1, ini.val=ini.val, 
                 unit="cm", par.list=FALSE, 
                 stand.fig=FALSE, angle=NULL, fig.opt=FALSE, 
                 control=list(reltol=1e-30, maxit=50000), 
                 np=2000 ) 

  n.loop <- 12
  Show   <- FALSE
  for(i in 1:n.loop){
    ini.val <- res0$par
    if(i==n.loop) Show <- TRUE
    print(paste(i, "/", n.loop, sep=""))
    res0 <- fitETE( x=x1, y=y1, ini.val=ini.val, 
                   unit="cm", par.list=FALSE, 
                   stand.fig=Show, angle=pi/4, fig.opt=Show,  
                   control=list(reltol=1e-30, maxit=50000), 
                   np=2000 )    
  }

  # The numerical values of the location and model parameters
  res0$par

  # The root-mean-square error (RMSE) between 
  #   the observed and predicted y values
  sqrt(res0$RSS/res0$sample.size)

  sqrt(sum((res0$y.stand.obs-res0$y.stand.pred)^2)/length(res0$y.stand.obs))

  # To calculate the volume of the egg
  VolumeETE(P=res0$par[4:7])

  # To calculate the surface area of the egg
  SurfaceAreaETE(P=res0$par[4:7])
}

graphics.off()
}