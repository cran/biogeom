\name{curveETE}
\encoding{UTF-8}
\alias{curveETE}
\title{
Drawing the Troscianko Curve Produced by the Explicit Troscianko Equation
}
\description{
\code{curveETE} is used to draw the Troscianko curve that is produced by the explicit Troscianko equation.
}
\usage{
curveETE(P, np = 5000, fig.opt = FALSE, deform.fun = NULL, Par = NULL,
        xlim = NULL, ylim = NULL, unit = NULL, main="")
}
\arguments{
  \item{P}{the three location parameters and the parameters of the explicit 
Troscianko equation.}
  \item{np}{the number of data points on the Troscianko curve.}
  \item{fig.opt}{an optional argument to draw the Troscianko curve.}
  \item{deform.fun}{the deformation function used to describe the deviation from a theoretical Troscianko curve.}
  \item{Par}{the parameter(s) of the deformation function.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Troscianko curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Troscianko curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Troscianko curve.}
  \item{main}{the main title of the figure.}
}
\details{
The first three elements of \code{P} are location parameters. The first two are the planar coordinates of the transferred origin, 
and the third is the angle between the major axis of the curve and the \eqn{x}-axis. Here, the major axis is a straight line through 
the two ends of an egg's profile (i.e., the mid-line of the egg's profile). The other arguments in \code{P} 
(except these first three location parameters) should correspond to those of \code{P} in \code{ETE}. 
\code{deform.fun} should take the form as: \code{deform.fun <- function(Par, z){\dots}}, where \code{z} is 
a two-dimensional matrix related to the \eqn{x} and \eqn{y} values. 
And the return value of \code{deform.fun} should be a \code{list} with two variables \code{x} and \code{y}.
}
\value{
  \item{x}{the \eqn{x} coordinates of the Troscianko curve.}
  \item{y}{the \eqn{y} coordinates of the Troscianko curve.}
}
\note{
When the rotation angle is zero (i.e., the third element in \code{P} is zero), \code{np} data points 
are distributed counterclockwise on the Troscianko curve from the rightmost end of the egg's profile to itself. 
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

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Wang, L., Quinn, B.K., Gielis, J. (2023) A new program to estimate the parameters 
of Preston's equation, a general formula for describing the egg shape of birds. \emph{Symmetry} 
15, 231. \doi{10.3390/sym15010231}
}
\seealso{
\code{\link{ETE}}, \code{\link{fitETE}}
}
\examples{
  Para1 <- c(0, 0, 0, 2.25, -0.377, -0.29, -0.16)
  curveETE(P=Para1, fig.opt=TRUE)

  # There is an example that introduces a deformation function in the egg-shape equation
  myfun <- function(Par, z){
    x  <- z[,1]
    y  <- z[,2]
    k1 <- Par[1]
    k2 <- Par[2]
    y  <- y - k1*(y+k2)^2
    list(x=x, y=y)
  }
  deform.op <- curveETE(P=Para1, np=5000, fig.opt=TRUE, deform.fun=myfun, Par=c(0.05, 8))

  graphics.off()
}