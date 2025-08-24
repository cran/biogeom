\name{fitAM}
\encoding{UTF-8}
\alias{fitAM}
\title{
Fitting the 2-D Profiles of Apical Meristems
}
\description{
\code{fitAM} is used to estimate the parameters of any non-linear equation 
the user defines for describing the 2-D profiles of shoot or root apical meristems.
}
\usage{
fitAM(expr, x, y, ini.val, method = "Nelder-Mead", control = list(),  
      lower = -Inf, upper = Inf, par.list = FALSE, stand.fig = TRUE, fig.opt = FALSE, 
      angle = NULL, xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{expr}{the non-linear equation for describing the 2-D profiles of shoot or root apical meristems.}
  \item{x}{the \eqn{x} coordinates of an apical meristem (AM)'s profile.}
  \item{y}{the \eqn{y} coordinates of an AM's profile.}
  \item{ini.val}{the list of initial values for the model parameters.}
  \item{method}{an optional argument to select an optimization method.}
  \item{control}{the list of control parameters for using the \code{\link[stats]{optim}} function in package \pkg{stats}.}
  \item{lower}{the lower bounds on the variables for the \code{L-BFGS-B} algorithm.}
  \item{upper}{the upper bounds on the variables for the \code{L-BFGS-B} algorithm.}
  \item{par.list}{an optional argument to show the list of parameters on the screen.}
  \item{stand.fig}{an optional argument to draw the observed and predicted profiles of an AM at the standard state 
(i.e., the origin for indicating the curve predicted by \code{expr} is located at (0, 0), 
and the \eqn{x}-axis of the predicted curve at this state is aligned with the \eqn{x}-axis).}
  \item{fig.opt}{an optional argument to draw the observed and predicted profiles of an AM 
at arbitrary angle between the rotated \eqn{x}-axis and the actual \eqn{x}-axis.}
  \item{angle}{the angle between the rotated \eqn{x}-axis and the actual \eqn{x}-axis, which can be defined by the user.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the predicted AM's profile curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the predicted AM's profile curve.}
  \item{unit}{the unit of the \eqn{x}-axis and the \eqn{y}-axis when showing the predicted AM's profile curve.}
  \item{main}{the main title of the figure.}
}
\details{
Here, the rotated \eqn{x}-axis means that for a scanned (photographed) AM the \eqn{x}-axis of the AM's profile  
curve at the standard state predicted by \code{expr} might deviate from the actual horizontal axis. 
The \code{Nelder-Mead} algorithm (Nelder and Mead, 1965) and the optimization method 
(referred to as \code{L-BFGS-B}) proposed by Byrd et al. (1995) in which each variable can be given 
a lower and/or upper bound can be selected to carry out the optimization of minimizing 
the residual sum of squares (RSS) between the observed and predicted \eqn{y} values. 
The \code{\link[stats]{optim}} function in package \pkg{stats} was used to carry out 
the \code{Nelder-Mead} algorithm and the \code{L-BFGS-B} algorithm. 
When \code{angle = NULL}, the observed AM's profile curve will be shown at its initial angle in the scanned image; 
when \code{angle} is a numerical value (e.g., \eqn{\pi/4}) defined by the user, 
it indicates that the \eqn{x}-axis at the standard state 
is rotated by the amount (\eqn{\pi/4}) counterclockwise from the actutal \eqn{x}-axis.
}
\value{
  \item{par}{the estimates of the model parameters.}
  \item{r.sq}{the coefficient of determination between the observed and predicted \eqn{y} values 
on the AM's profile curve.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values 
on the AM's profile curve.}
  \item{x.stand.obs}{the observed \eqn{x} coordinates of the points 
on the AM's profile curve at the standard state.}
  \item{x.stand.pred}{the predicted \eqn{x} coordinates of the points 
on the AM's profile curve at the standard state.} 
  \item{y.stand.obs}{the observed \eqn{y} coordinates of the points 
on the AM's profile curve at the standard state.}
  \item{y.stand.pred}{the predicted \eqn{y} coordinates of the points 
on the AM's profile curve at the standard state.} 
  \item{x.obs}{the observed \eqn{x} coordinates of the points on the AM's profile curve 
at the transferred polar angles as defined by the user.}
  \item{x.pred}{the predicted \eqn{x} coordinates of the points on the AM's profile curve 
at the transferred polar angles as defined by the user.}   
  \item{y.obs}{the observed \eqn{y} coordinates of the points on the AM's profile curve 
at the transferred polar angles as defined by the user.}
  \item{y.pred}{the predicted \eqn{y} coordinates of the points on the AM's profile curve 
at the transferred polar angles as defined by the user.}       
}
\note{
In the arguments, \code{expr} can be flexibly defined by the user, but it requires taking the form as 
\code{myfun <- function(P, x){...}}, where \code{P} represents the parameter vector, and \code{x} 
is the independent variable. In the outputs, the data length of the predicted values is the same as that of 
the observations. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Byrd, R.H., Lu, P., Nocedal, J., Zhu, C. (1995) A limited memory algorithm for bound constrained 
optimization. \emph{SIAM Journal on Scientific Computing} 16, 1190\eqn{-}1208. \doi{10.1137/0916069}

Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2023) 
A simple way to calculate the volume and surface area of avian eggs. 
\emph{Annals of the New York Academy of Sciences} 1524, 118\eqn{-}131. \doi{10.1111/nyas.15000}
}
\seealso{
\code{\link{PlanCoor}}, \code{\link{SAMs}}, \code{\link{SurfaceAreaAM}}, \code{\link{VolumeAM}}
}
\examples{

#### See Shi et al. (2025) for details #########################################
# Shi, P., Liu, X., Gielis, J., Beirinckx, B., Niklas, K.J. (2025) 
#     Comparison of six non-linear equations in describing the 2-D 
#     profiles of apical meristems. American Journal of Botany (under review).
################################################################################

data(SAMs)

uni.sam <- sort( unique(SAMs$Genus) )
ind     <- 2
Data    <- SAMs[SAMs$Genus==uni.sam[ind], ]
x0      <- Data$x
y0      <- Data$y

Res1    <- adjdata(x0, y0, ub.np=200, times=1.2, len.pro=1/20)
X       <- Res1$x
Y       <- Res1$y

dev.new()
plot( X, Y, pch=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(paste(italic(x), " (unitless)", sep="")), 
      ylab=expression(paste(italic(y), " (unitless)", sep="")) )

x1 <- X
y1 <- Y

if(TRUE){

  # The code results in an upward-opening profile curve with 
  #   all y-values less than zero, facilitating data fitting by 
  #   providing suitable initial values for model parameters.

  y1 <- min(Y)-Y
  x1 <- X-mean(X)  
   
  # To normalize the profile coordinates, 
  #   and make x-coordinates range between -1 and 1  
  x2 <- x1/max(abs(x1))
  y2 <- y1/max(abs(x1))

  rm(x1)
  rm(y1)

  x1 <- x2
  y1 <- y2
}


#### Fitting the catenary equation ###############################################
myfun1 <- function(P, x){
  a <- P[1]
  return( a * cosh(x/a) )
}

x0.ini    <- 0
y0.ini    <- -min(y1)
theta.ini <- 0 
a.ini     <- c( 0.75*abs(min(y1)), abs(min(y1)), 1.25*abs(min(y1)) )
ini.val1  <- list(x0.ini, y0.ini, theta.ini, a.ini )

Res1 <- fitAM( myfun1, x=x1, y=y1, ini.val=ini.val1, 
               control=list(reltol=1e-20, maxit=20000), 
               par.list=FALSE, stand.fig=FALSE, fig.opt=TRUE, angle=NULL, 
               unit="unitless", main="Catenary equation" )
###################################################################################


#### Fitting the parabolic equation ###############################################
myfun2 <- function(P, x){
  alpha1 <- P[1]
  alpha2 <- P[2]
  alpha1*x + alpha2*x^2 
}

x0.ini    <- 0
y0.ini    <- -abs(min(y1))
theta.ini <- 0 
beta1.ini <- 1e-3
beta2.ini <- c(1, 5, 10, 15)
ini.val2  <- list(x0.ini, y0.ini, theta.ini, beta1.ini, beta2.ini)

Res2 <- fitAM( myfun2, x=x1, y=y1, ini.val=ini.val2, 
               control=list(factr=1e-12, maxit=20000), 
               method="L-BFGS-B", lower=c(-Inf, -Inf, -pi/2/4, -Inf, 0), 
               upper=c(Inf, Inf, pi/2/4, Inf, Inf), 
               par.list=FALSE, stand.fig=FALSE, fig.opt=TRUE, angle=NULL, 
               unit="unitless", main="Parabolic equation" )
###################################################################################


#### Fitting the hybrid catenary-parabolic equation ###############################
myfun3 <- function(P, x){
  alpha <- P[1]
  beta  <- P[2]
  gamma <- P[3]
  alpha * cosh(beta * x) + gamma * x^2 - alpha
}

x0.ini    <- 0
y0.ini    <- -abs(min(y1))
theta.ini <- 0 
alpha.ini <- c(0.1, 0.5, 1, 2.5, 5)
beta.ini  <- 1/abs(min(y1))
gamma.ini <- 1
ini.val3  <- list(x0.ini, y0.ini, theta.ini, alpha.ini, beta.ini, gamma.ini)

Res3 <- fitAM( myfun3, x=x1, y=y1, ini.val=ini.val3,  
               control=list(reltol=1e-30, maxit=50000),           
               par.list=FALSE, stand.fig=FALSE, fig.opt=TRUE, angle=NULL,               
               unit="unitless", main="Hybrid catenary-parabolic equation" )
###################################################################################


#### Fitting the performance equation #############################################
myfun4 <- function(P, x){
  c    <- P[1]
  K1   <- P[2]
  K2   <- P[3]
  xmin <- 0
  xmax <- P[4]
  # x[x > xmax] <- xmax
  # x[x < xmin] <- xmin
  -c * ( 1-exp(-K1*(x-xmin)) )*( 1-exp(K2*(x-xmax)) )
}

x0.ini    <- min(x1)
y0.ini    <- 0
theta.ini <- -pi/12 
c.ini     <- abs(min(y1))
K1.ini    <- 1
K2.ini    <- 1
xmax.ini  <- max(x1)*2
ini.val4  <- list(x0.ini, y0.ini, theta.ini, c.ini, K1.ini, K2.ini, xmax.ini)

Res4 <- fitAM( myfun4, x=x1, y=y1, ini.val=ini.val4, method="L-BFGS-B",
               lower=c(-Inf, -Inf, -pi/2/4, 0, 0, 0, 0), 
               upper=c(Inf, Inf, pi/2/4, Inf, Inf, Inf, Inf),   
               control=list(factr=1e-12, maxit=20000),           
               par.list=FALSE, stand.fig=FALSE, fig.opt=TRUE, angle=NULL,
               unit="unitless", main="Performance equation" )
###################################################################################


#### Fitting the superparabolic equation equation #################################
myfun5 <- function(P, x){
  beta1 <- P[1]
  beta2 <- P[2]
  beta1 * abs(x)^beta2
}

x0.ini    <- 0
y0.ini    <- -max( y1 )  
theta.ini <- 0
beta1.ini <- max(x1)/2
beta2.ini <- c(1.5, 2, 2.5)
ini.val5  <- list(x0.ini, y0.ini, theta.ini, beta1.ini, beta2.ini)

Res5 <- fitAM( myfun5, x=x1, y=y1, ini.val=ini.val5, 
               control=list(reltol=1e-20, maxit=20000), 
               par.list=FALSE, stand.fig=FALSE, fig.opt=TRUE, angle=NULL,
               unit="unitless", main="Superparabolic equation" )
###################################################################################


#### Fitting the superellipse equation equation ###################################
myfun6 <- function(P, x){
  A <- P[1]
  B <- P[2]
  n <- P[3]
  -B*(1-abs(x/A)^n)*(1/n)
}

x0.ini    <- 0
y0.ini    <- 0
theta.ini <- c(0, pi/2) 
A.ini     <- max(x1)
B.ini     <- -min(y1)
n.ini     <- c(1, 2, 3)
ini.val6  <- list(x0.ini, y0.ini, theta.ini, A.ini, B.ini, n.ini)
Res6 <- fitAM( myfun6, x=x1, y=y1, ini.val=ini.val6, 
               control=list(reltol=1e-20, maxit=20000), 
               par.list=FALSE, stand.fig = FALSE, fig.opt=TRUE, angle=NULL,
               unit="unitless", main="Superellipse equation" )
###################################################################################

graphics.off()
}