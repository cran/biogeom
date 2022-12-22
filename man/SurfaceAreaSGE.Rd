\name{SurfaceAreaSGE}
\encoding{UTF-8}
\alias{SurfaceAreaSGE}
\title{
Calculation of the Surface Area of An Egg Based on the Simplified Gielis Equation 
}
\description{
\code{SurfaceAreaSGE} is used to calculate the surface area of an egg that follows the simplified Gielis equation.
}
\usage{
SurfaceAreaSGE(P, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{P}{the parameters of the simplified Gielis equation, including \eqn{a}, \eqn{n_{1}}, and \eqn{n_{2}}.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
}
\details{
The formula of the surface area (\eqn{S}) of an egg based on the simplified Gielis equation is:
\deqn{S(\varphi)=2\,\pi\int_{0}^{\pi}\sin{\left(\varphi\right)}\ r\,\sqrt{r^2+\left(\frac{dr}{d\varphi}\right)^2}\,d\varphi,}

where the polar raidus (\eqn{r}) is the function of the polar angle (\eqn{\varphi}):

\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}\right)^{-\frac{1}{n_{1}}},} 

namely the simplified Gielis equation (i.e., \code{\link{GE}}) with arguments 
\code{simpver = 1} and \code{m = 1}.  
}
\note{
The argument \code{P} in the \code{\link{SurfaceAreaSGE}} function only has 
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
\code{\link{DSGE}}, \code{\link{fitGE}}, \code{\link{GE}}, \code{\link{VolumeSGE}}
}
\examples{
  Par7 <- c(1.124, 14.86, 49.43)
  SurfaceAreaSGE(P = Par7)
}