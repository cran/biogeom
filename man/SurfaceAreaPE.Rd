\name{SurfaceAreaPE}
\encoding{UTF-8}
\alias{SurfaceAreaPE}
\title{
Calculation of the Surface Area of An Egg Based on the Preston Equation 
}
\description{
\code{SurfaceAreaPE} is used to calculate the surface area of an egg that follows the Preston equation.
}
\usage{
SurfaceAreaPE(P, simpver = NULL, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{P}{the parameters of the Preston equation or one of its simplified versions.}
  \item{simpver}{an optional argument to use the simplified version of the explicit Preston equation.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
}
\details{
The formula of the surface area of an egg based on the explicit Preston equation or one of its simplified versions is:
\deqn{S(x)=2\,\pi\int_{-a}^{a}y\,\sqrt{1+\left(\frac{dy}{dx}\right)^2}\,dx}
Here, \eqn{y} denotes the explicit Preston equation (i.e., \code{\link{EPE}}).  
When \code{simpver = NULL}, \code{P} has five parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, \eqn{c_{2}}, and \eqn{c_{3}}; 
when \code{simpver = 1}, \code{P} has four parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, and \eqn{c_{2}}; 
when \code{simpver = 2}, \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{1}}; 
when \code{simpver = 3}, \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{2}}.  
}
\note{
The \code{\link{SurfaceAreaPE}} function has the same parameters, \code{P}, as those in the 
\code{\link{PE}} function. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
A manuscript for \emph{Annals of the New York Academy of Sciences} is 
currently under preparation by Shi et al.
}
\seealso{
\code{\link{DEPE}}, \code{\link{EPE}}, \code{\link{PE}}, \code{\link{VolumePE}}
}
\examples{
  Par3 <- c(4.27, 2.90, 0.0868, 0.0224, -0.0287)
  SurfaceAreaPE(P = Par3, simpver = NULL)
}