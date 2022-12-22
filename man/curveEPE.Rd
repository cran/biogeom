\name{curveEPE}
\encoding{UTF-8}
\alias{curveEPE}
\title{
Drawing the Preston Curve Produced by the the Explicit Preston Equation
}
\description{
\code{curveEPE} is used to draw the Preston curve that is produced by the explicit Preston equation.
}
\usage{
curveEPE(P, np = 5000, simpver = NULL,
        fig.opt = FALSE, deform.fun = NULL, Par = NULL,
        xlim = NULL, ylim = NULL, unit = NULL, main="")
}
\arguments{
  \item{P}{the three location parameters and the parameters of the explicit 
Preston equation or one of its simplified versions.}
  \item{np}{the number of data points on the Preston curve.}
  \item{simpver}{an optional argument to use the simplfied version of the explicit Preston equation.}
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
the two ends of an egg's profile (i.e., the mid-line of the egg's profile). The other arguments in \code{P} 
(except these first three location parameters), and \code{simpver} should correspond to those of \code{P} in \code{EPE}. 
\code{deform.fun} should take the form as: \code{deform.fun <- function(Par, z){\dots}}, where \code{z} is 
a two-dimensional matrix related to the \eqn{x} and \eqn{y} values. 
And the return value of \code{deform.fun} should be a \code{list} with two variables \code{x} and \code{y}.
}
\value{
  \item{x}{the \eqn{x} coordinates of the Preston curve.}
  \item{y}{the \eqn{y} coordinates of the Preston curve.}
}
\note{
When the rotation angle is zero (i.e., the third element in \code{P} is zero), \code{np} data points 
are distributed counterclockwise on the Preston curve from the rightmost end of the egg's profile to itself. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2022) 
Can we simply calculate the volume and surface area of an egg? Evidence from six species of poultry. 
\emph{Annals of the New York Academy of Sciences} Under preparation.

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Wang, L., Quinn, B.K., Gielis, J. (2022) A new program to estimate the parameters 
of Preston's equation, a general formula for describing the egg shape of birds. \emph{Symmetry} 
Under review. 

Todd, P.H., Smart, I.H.M. (1984) The shape of birds' eggs. \emph{Journal of Theoretical Biology} 
106, 239\eqn{-}243. \doi{10.1016/0022-5193(84)90021-3}
}
\seealso{
\code{\link{EPE}}, \code{\link{fitEPE}}, \code{\link{lmPE}}, \code{\link{PE}}, \code{\link{TSE}}
}
\examples{
  Para1 <- c(0, 0, 0, 10, 6, 0.325, -0.0415)
  curveEPE(P=Para1, simpver=1, fig.opt=TRUE)
  Para2 <- c(0, 0, pi, 10, 6, -0.325, -0.0415)
  curveEPE(P=Para2, simpver=1, fig.opt=TRUE)

  Para3 <- c(0, 0, 0, 10, 6, 0.325, -0.0415, 0.2)
  curveEPE(P=Para3, simpver=NULL, fig.opt=TRUE)
  Para4 <- c(0, 0, pi, 10, 6, -0.325, -0.0415, 0.2)
  curveEPE(P=Para4, simpver=NULL, fig.opt=TRUE)

  Para5 <- c(0, 0, pi/4, 10, 6, 0.325, -0.0415)
  curveEPE(P=Para5, simpver=1, 
          fig.opt=TRUE, main="A rotated egg shape")

  # There is an example that introduces a deformation function in the egg-shape equation
  myfun <- function(Par, z){
    x  <- z[,1]
    y  <- z[,2]
    k1 <- Par[1]
    k2 <- Par[2]
    y  <- y - k1*(y+k2)^2
    list(x=x, y=y)
  }
  deform.op <- curveEPE(P=Para1, np=5000, simpver=1, 
                       fig.opt=TRUE, deform.fun=myfun, Par=c(0.05, 8))


  graphics.off()
}