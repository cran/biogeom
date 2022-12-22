\name{DSGE}
\encoding{UTF-8}
\alias{DSGE}
\title{
Calculation of the First-Order Derivative of the Simplified Gielis Equation  
}
\description{
\code{DSGE} is used to calculate the first-order derivative of the simplified Gielis equation at a given \eqn{\varphi}-value.
}
\usage{
DSGE(P, phi)
}
\arguments{
  \item{P}{the parameters of the simplified Gielis equation, including \eqn{a}, \eqn{n_{1}}, and \eqn{n_{2}}.}
  \item{phi}{the \eqn{\varphi}-value used in the simplified Gielis equation.}
}
\details{
The first-order derivative of the simplified Gielis equation with arguments \code{simpver = 1} and \code{m = 1} at a given \eqn{\varphi}-value is: 
\deqn{g(x) = \frac{a}{4}\,\frac{n_{2}}{n_{1}}\,\left[\left(\cos{\frac{\varphi}{4}}\right)^{n_{2}-1}\left(\sin{\frac{\varphi}{4}}\right)-\left(\sin{\frac{\varphi}{4}}\right)^{n_{2}-1}\left(\cos{\frac{\varphi}{4}}\right)\right]\ \left[\left(\cos{\frac{\varphi}{4}}\right)^{n_{2}}+\left(\sin{\frac{\varphi}{4}}\right)^{n_{2}}\right]^{-\frac{1}{n_{1}}-1}, }

where \code{P} has three parameters: \eqn{a}, \eqn{n_{1}}, and \eqn{n_{2}}.
}
\note{
The argument \code{P} in the \code{\link{DSGE}} function only has 
the three parameters: \eqn{a}, \eqn{n_{1}}, and \eqn{n_{2}}. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Chen, Z. (2012) Volume and area of revolution under polar coordinate system. 
\emph{Studies in College Mathematics} 15(6), 9\eqn{-}11.

Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2022) 
Can we simply calculate the volume and surface area of an egg? Evidence from six species of poultry. 
\emph{Annals of the New York Academy of Sciences} Under preparation.
}
\seealso{
\code{\link{GE}}, \code{\link{fitGE}}, \code{\link{SurfaceAreaSGE}}
}
\examples{
  Par7 <- c(1.124, 14.86, 49.43)
  phi1 <- seq(0, pi, len=2000)
  g1   <- DSGE(P=Par7, phi=phi1)

  dev.new()
  plot(phi1, g1, type="l", col=4, cex.lab=1.5, cex.axis=1.5,
       xlab=expression(italic(varphi)), 
       ylab=expression(paste(italic(g), "(", italic(varphi), ")", sep="")))

  graphics.off()
}