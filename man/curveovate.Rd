\name{curveovate}
\encoding{UTF-8}
\alias{curveovate}
\title{
Drawing the Ovate Leaf-Shape Curve
}
\description{
\code{curveovate} is used to draw the ovate leaf-shape curve.
}
\usage{
curveovate(expr, P, x, fig.opt = FALSE, 
           deform.fun = NULL, Par = NULL,
           xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{
  \item{expr}{the simplified version 1 of a performance equation.}
  \item{P}{the three location parameters and the parameters of the simplified version 1 of a performance equation.}
  \item{x}{the given \eqn{x} values to draw the ovate leaf-shape curve.}
  \item{fig.opt}{an optional argument to draw the ovate leaf-shape curve.}
  \item{deform.fun}{the deformation function used to describe the deviation 
from a theoretical ovate leaf-shape curve.}
  \item{Par}{the parameter(s) of the deformation function.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the ovate leaf-shape curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the ovate leaf-shape curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the ovate leaf-shape curve.}
  \item{main}{the main title of the figure.}
}
\details{
\code{P} has seven elements: three location parameters, and four model parameters, 
i.e., \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, 
\eqn{x_{\mathrm{max}}}, and \eqn{\delta}. This means that \code{expr} is limited to 
be the simplified version 1 (where \eqn{x_{\mathrm{min}} = 0})  
in \code{\link{MbetaE}}, \code{\link{MBriereE}}, and \code{\link{MLRFE}}. 
The first three elements of \code{P} are location parameters, among which the first two 
are the planar coordinates of the transferred origin, 
and the third is the angle between the major axis of the curve and the \eqn{x}-axis. 
\code{deform.fun} should take the form as: \code{deform.fun <- function(Par, z){\dots}}, where \code{z} is 
a two-dimensional matrix related to the \eqn{x} and \eqn{y} values. 
And the return value of \code{deform.fun} should be a \code{list} with two variables \code{x} and \code{y}.
}
\value{
  \item{x}{the \eqn{x} coordinates of the ovate leaf-shape curve.}
  \item{y}{the \eqn{y} coordinates of the ovate leaf-shape curve.} 
}
\note{
The number of elements in \code{P} here has additional three location parameters than that 
in \code{\link{MbetaE}}, \code{\link{MBriereE}}, and \code{\link{MLRFE}}.
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
\code{\link{areaovate}}, \code{\link{fitovate}}, \code{\link{MbetaE}}, \code{\link{MBriereE}}, \code{\link{MLRFE}}
}
\examples{
P1  <- c(1, 1, pi/4, 2, 3, 10, 4)
RE1 <- curveovate(MLRFE, P=P1, x=seq(0, 10, by=0.1), fig.opt=TRUE)
RE2 <- curveovate(MbetaE, P=P1, x=seq(0, 10, by=0.1), fig.opt=TRUE)

dev.new()
plot(RE1$x, RE1$y, cex.lab=1.5, cex.axis=1.5, type="l", 
  xlab=expression(italic("x")), ylab=expression(italic("y"))) 
lines(RE2$x, RE2$y, col=4)
}