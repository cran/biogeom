\name{curvePE}
\encoding{UTF-8}
\alias{curvePE}
\title{
Drawing the Preston Curve
}
\description{
\code{curvePE} is used to draw the Preston curve.
}
\usage{
curvePE(P, zeta = seq(0, 2*pi, len = 2000), simpver = NULL,
        fig.opt = FALSE, deform.fun = NULL, Par = NULL,
        xlim = NULL, ylim = NULL, unit = NULL, main="")
}
\arguments{
  \item{P}{the three location parameters and the parameters of the original 
Preston equation or one of its simplified versions.}
  \item{zeta}{the given angles at which we want to draw the Preston curve.}
  \item{simpver}{an optional argument to use the simplfied version of the original Preston equation.}
  \item{fig.opt}{an optional argument to draw the Preston curve.}
  \item{deform.fun}{the deformation function used to describe the deviation from a theoretical Preston curve.}
  \item{Par}{the parameter(s) of the deformation function.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Preston curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Preston curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Preston curve.}
  \item{main}{the main title of the figure.}
}
\details{
The first three elements of \code{P} are location parameters. The first two are the planar coordinates of the transferred origin, 
and the third is the angle between the major axis of the curve and the \eqn{x}-axis. Here, the major axis is a straight line through 
the midpoint of the egg length, and is perpendicular to the egg length axis. The other arguments in \code{P} 
(except these first three location parameters), and \code{simpver} should correspond to \code{PE}. 
\code{deform.fun} should take the form as: \code{deform.fun <- function(Par, z){\dots}}, where \code{z} is 
a two-dimensional matrix related to the \eqn{x} and \eqn{y} values. 
And the return value of \code{deform.fun} should be a \code{list} with two variables \code{x} and \code{y}.
}
\value{
  \item{x}{the \eqn{x} coordinates of the Preston curve corresponding to the given angles \code{zeta}.}
  \item{y}{the \eqn{y} coordinates of the Gielis curve corresponding to the given angles \code{zeta}.}
  \item{r}{the distances of the Preston curve from the origin corresponding to the given angles \code{zeta}.}
}
\note{
\eqn{\zeta} is NOT the polar angle corresponding to \eqn{r}, i.e., 
\deqn{y \neq r\,\mathrm{sin}\,\zeta} 
\deqn{x \neq r\,\mathrm{cos}\,\zeta} 
Let \eqn{\varphi} be the polar angle corresponding to \eqn{r}. We have: 
\deqn{ \zeta = \mathrm{arc\,sin}\frac{ r\ \mathrm{sin}\,\varphi }{a} }
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

Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14862}

Todd, P.H., Smart, I.H.M. (1984) The shape of birds' eggs. \emph{Journal of Theoretical Biology} 
106, 239\eqn{-}243. \doi{10.1016/0022-5193(84)90021-3}
}
\seealso{
\code{\link{fitPE}}, \code{\link{PE}}, \code{\link{lmPE}}, \code{\link{TSE}}
}
\examples{
  Para1 <- c(0, 0, 0, 10, 6, 0.325, -0.0415)
  curvePE(P=Para1, zeta=seq(0, 2*pi, len=2000), simpver=1, fig.opt=TRUE)
  Para2 <- c(0, 0, pi, 10, 6, -0.325, -0.0415)
  curvePE(P=Para2, zeta=seq(0, 2*pi, len=2000), simpver=1, fig.opt=TRUE)

  Para3 <- c(0, 0, 0, 10, 6, 0.325, -0.0415, 0.2)
  curvePE(P=Para3, zeta=seq(0, 2*pi, len=2000), simpver=NULL, fig.opt=TRUE)
  Para4 <- c(0, 0, pi, 10, 6, -0.325, -0.0415, 0.2)
  curvePE(P=Para4, zeta=seq(0, 2*pi, len=2000), simpver=NULL, fig.opt=TRUE)

  Para5 <- c(0, 0, pi/4, 10, 6, 0.325, -0.0415)
  curvePE(P=Para5, zeta=seq(0, 2*pi, len=2000), simpver=1, 
          fig.opt=TRUE, main="A rotated egg shape")

  # There is an example that introduces a deformation function in the egg-shape equation
  myfun <- function(Par, z){
    x  <- z[,1]
    y  <- z[,2]
    k1 <- Par[1]
    k2 <- Par[2]
    x  <- x - k1*(x+k2)^2
    list(x=x, y=y)
  }
  deform.op <- curvePE(P=Para1, zeta=seq(0, 2*pi, len=2000), simpver=1, 
                       fig.opt=TRUE, deform.fun=myfun, Par=c(0.05, 8))

}