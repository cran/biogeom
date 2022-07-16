\name{curveNRGE}
\encoding{UTF-8}
\alias{curveNRGE}
\title{
Drawing the Egg Shape Predicted by the Narushin-Romanov-Griffin Equation
}
\description{
\code{curveNRGE} is used to draw the egg shape predicted by the Narushin-Romanov-Griffin equation.
}
\usage{
curveNRGE(P, x, fig.opt = FALSE, deform.fun = NULL, 
          Par = NULL, xlim = NULL, ylim = NULL, unit = NULL, main = NULL)
}
\arguments{ 
  \item{P}{the three location parameters and the four parameters of the Narushin-Romanov-Griffin equation (Narushin et al., 2021).}
  \item{x}{the given \eqn{x} coordinates at which we want to draw the Narushin-Romanov-Griffin curve.}
  \item{fig.opt}{an optional argument to draw the Narushin-Romanov-Griffin curve.}
  \item{deform.fun}{the deformation function used to describe the deviation from a theoretical Narushin-Romanov-Griffin curve.}
  \item{Par}{the parameter(s) of the deformation function.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot the Narushin-Romanov-Griffin curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot the Narushin-Romanov-Griffin curve.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the Narushin-Romanov-Griffin curve.}
  \item{main}{the main title of the figure.}
}
\details{
The first three elements of \code{P} are location parameters. The first two are the planar coordinates of the transferred origin, 
and the third is the angle between the major axis of the curve and the \eqn{x}-axis. The other arguments in \code{P} 
should be the same as those in \code{\link{NRGE}}.  
\code{deform.fun} should take the form as: \code{deform.fun <- function(Par, z){\dots}}, where \code{z} is 
a two-dimensional matrix related to the \eqn{x} and \eqn{y} values. 
And the return value of \code{deform.fun} should be a \code{list} with two variables \code{x} and \code{y}.
}
\value{
  \item{x}{the \eqn{x} coordinates of the Narushin-Romanov-Griffin curve.}
  \item{y}{the \eqn{y} coordinates of the Narushin-Romanov-Griffin curve.}
}
\note{
The \eqn{x} coordinates of the Narushin-Romanov-Griffin curve are different from the given \eqn{x} coordinates 
(i.e., \code{x} in Arguments). The latter are twice the former. The \eqn{x} coordinates for the 
Narushin-Romanov-Griffin curve range from the maximum to the minimum on the \eqn{x}-axis (corresponding to positive \eqn{y} coordinates), 
and then range from the minimum to the maximum on the \eqn{x}-axis (corresponding to negative \eqn{y} coordinates). That is to say, 
the coordinates of \eqn{(x, y)} need to undergo a counterclockwise rotation to form the Narushin-Romanov-Griffin curve.       
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Narushin, V.G., Romanov, M.N., Griffin, D.K. (2021) Egg and math: introducing a universal formula for egg shape. 
\emph{Annals of the New York Academy of Sciences} 1505, 169\eqn{-}177. \doi{10.1111/nyas.14680}

Shi, P., Gielis, J., Niklas, K.J. (2022) Comparison of a universal (but complex) model for avian egg 
shape with a simpler model. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14799}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14862}

Tian, F., Wang, Y., Sandhu, H.S., Gielis, J., Shi, P. (2020) Comparison of seed morphology of two ginkgo cultivars. 
\emph{Journal of Forestry Research} 31, 751\eqn{-}758. \doi{10.1007/s11676-018-0770-y}
}
\seealso{
\code{\link{fitNRGE}}, \code{\link{NRGE}}
}
\examples{
PA   <- c(1, 1, pi/4, 11.5, 7.8, 1.1, 5.6)
xA   <- seq(-11.5/2, 11.5/2, len=2000)
resA <- curveNRGE(PA, xA, fig.opt=TRUE)
xB   <- seq(-11.5/2, 11.5/2, len=100)
resB <- curveNRGE(PA, xB, fig.opt=TRUE, xlim=c(-6, 6), 
                  ylim=c(-6, 6), main="A pear-shaped egg") 
cbind(resB$x, resB$y)
}




