\name{curveGE}
\encoding{UTF-8}
\alias{curveGE}
\title{
Drawing the Gielis Curve
}
\description{
\code{curveGE} is used to draw the Gielis curve.
}
\usage{
curveGE(expr, P, phi = seq(0, 2*pi, len = 2000),
        m = 1, simpver = NULL, nval = 1,
        fig.opt = FALSE, deform.fun = NULL, Par = NULL,
        xlim = NULL, ylim = NULL, unit = NULL, main="")
}
\arguments{
  \item{expr}{the original (or twin) Gielis equation or one of its simplified versions.}
  \item{P}{the three location parameters and the parameters of the original (or twin) 
Gielis equation or one of its simplified versions.}
  \item{phi}{the given polar angles at which we want to draw the Gielis curve.}
  \item{m}{the given \eqn{m} value that determines the number of angles of the Gielis curve within \eqn{[0, 2\pi)}.}
  \item{simpver}{an optional argument to use the simplfied version of the original (or twin) Gielis equation.}
  \item{nval}{the specified value for \eqn{n_{1}} or \eqn{n_{2}} or \eqn{n_{3}} in the simplified versions.}
  \item{fig.opt}{an optional argument to draw the Gielis curve.}
  \item{deform.fun}{the deformation function used to describe the deviation from a theoretical Gielis curve.}
  \item{Par}{the parameter(s) of the deformation function.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Gielis curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Gielis curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Gielis curve.}
  \item{main}{the main title of the figure.}
}
\details{
The first three elements of \code{P} are location parameters. The first two are the planar coordinates of the transferred polar point, 
and the third is the angle between the major axis of the curve and the \eqn{x}-axis. The other arguments in \code{P} 
(except these first three location parameters), \code{m}, \code{simpver}, and \code{nval} should correspond 
to \code{expr} (i.e., \code{GE} or \code{TGE}).  
Please note the differences in the simplified version number and the 
number of parameters between \code{GE} and \code{TGE}. 
\code{deform.fun} should take the form as: \code{deform.fun <- function(Par, z){\dots}}, where \code{z} is 
a two-dimensional matrix related to the \eqn{x} and \eqn{y} values. 
And the return value of \code{deform.fun} should be a \code{list} with two variables \code{x} and \code{y}.
}
\value{
  \item{x}{the \eqn{x} coordinates of the Gielis curve corresponding to the given polar angles \code{phi}.}
  \item{y}{the \eqn{y} coordinates of the Gielis curve corresponding to the given polar angles \code{phi}.}
  \item{r}{the polar radii of the Gielis curve corresponding to the given polar angles \code{phi}.}
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
and abstract shapes. \emph{American Journal of Botany} 90, 333-338. \doi{10.3732/ajb.90.3.333}

Li, Y., Quinn, B.K., Gielis, J., Li, Y., Shi, P. (2022) Evidence that supertriangles exist 
in nature from the vertical projections of \emph{Koelreuteria paniculata} fruit. \emph{Symmetry} 
14, 23. \doi{10.3390/sym14010023}

Shi, P., Gielis, J., Niklas, K.J. (2022) Comparison of a universal (but complex) model for avian egg 
shape with a simpler model. \emph{Annals of the New York Academy of Sciences} 1514, 34\eqn{-}42. \doi{10.1111/nyas.14799}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Ratkowsky, D.A., Gielis, J. (2020) The generalized Gielis geometric
equation and its application. \emph{Symmetry} 12, 645. \doi{10.3390/sym12040645}

Shi, P., Xu, Q., Sandhu, H.S., Gielis, J., Ding, Y., Li, H., Dong, X. (2015) Comparison
of dwarf bamboos (\emph{Indocalamus} sp.) leaf parameters to determine relationship
between spatial density of plants and total leaf area per plant. \emph{Ecology and Evolution} 5,
4578-4589. \doi{10.1002/ece3.1728}
}
\seealso{
\code{\link{areaGE}}, \code{\link{fitGE}}, \code{\link{GE}}, \code{\link{TGE}}
}
\examples{
GE.par  <- c(2, 1, 4, 6, 3)
phi.vec <- seq(0, 2*pi, len=2000)
r.theor <- GE(P=GE.par, phi=phi.vec, m=5)

dev.new()
plot( phi.vec, r.theor, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic(phi)), ylab=expression(italic("r")),
      type="l", col=4 ) 

curve.par <- c(1, 1, pi/4, GE.par)
GE.res    <- curveGE(GE, P=curve.par, fig.opt=TRUE, deform.fun=NULL, Par=NULL, m=5)
# GE.res$r

GE.res    <- curveGE( GE, P=c(0, 0, 0, 2, 4, 20), m=1, simpver=1, fig.opt=TRUE )
# GE.res$r

GE.res    <- curveGE( GE, P=c(1, 1, pi/4, 2, 1, 3), m=5, simpver=1, fig.opt=TRUE )
# GE.res$r

GE.res    <- curveGE( GE, P=c(1, 1, pi/4, 2, 1, 3), m=2, simpver=1, fig.opt=TRUE )
# GE.res$r

GE.res    <- curveGE( GE, P=c(1, 1, pi/4, 2, 0.05), m=1, simpver=2, fig.opt=TRUE )
# GE.res$r

GE.res    <- curveGE( GE, P=c(1, 1, pi/4, 2), m=4, simpver=3, nval=2, fig.opt=TRUE )
# GE.res$r

GE.res    <- curveGE( GE, P=c(1, 1, pi/4, 2, 0.6), m=4, simpver=8, nval=2, fig.opt=TRUE )
# GE.res$r

graphics.off()
}