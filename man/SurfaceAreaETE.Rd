\name{SurfaceAreaETE}
\encoding{UTF-8}
\alias{SurfaceAreaETE}
\title{
Calculation of the Surface Area of An Egg Based on the Explicit Troscianko Equation 
}
\description{
\code{SurfaceAreaETE} is used to calculate the surface area of an egg that follows the explicit Troscianko equation.
}
\usage{
SurfaceAreaETE(P, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{P}{the parameters of the explicit Troscianko equation.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
}
\details{
The formula of the surface area (\eqn{S}) of an egg based on the explicit Troscianko equation is:
\deqn{S(x)=2\,\pi\int_{-a}^{a}y\,\sqrt{1+\left(\frac{dy}{dx}\right)^2}\,dx,}
where \eqn{y} denotes the explicit Troscianko equation (i.e., \code{\link{ETE}}), and 
\eqn{a} denotes half the egg's length.  
}
\note{
The argument \code{P} in the \code{\link{SurfaceAreaETE}} function has the same parameters, as those in the 
\code{\link{ETE}} function. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Narushin, V.G., Romanov, M.N., Mishra, B., Griffin, D.K. (2022) Mathematical progression of 
avian egg shape with associated area and volume determinations. 
\emph{Annals of the New York Academy of Sciences} 1513, 65\eqn{-}78. \doi{10.1111/nyas.14771}

Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2023) 
A simple way to calculate the volume and surface area of avian eggs. 
\emph{Annals of the New York Academy of Sciences} 1524, 118\eqn{-}131. \doi{10.1111/nyas.15000}
}
\seealso{
\code{\link{DETE}}, \code{\link{ETE}}, \code{\link{fitETE}}, \code{\link{VolumeETE}}
}
\examples{
  Par5 <- c(2.25, -0.38, -0.29, -0.16)
  SurfaceAreaETE(P = Par5)
}