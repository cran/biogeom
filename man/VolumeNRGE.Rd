\name{VolumeNRGE}
\encoding{UTF-8}
\alias{VolumeNRGE}
\title{
Calculation of the Volume of An Egg Based on the Narushin-Romanov-Griffin Equation
}
\description{
\code{VolumeNRGE} is used to calculate the volume of an egg that follows the Narushin-Romanov-Griffin equation.
}
\usage{
VolumeNRGE(P, subdivisions = 100L,
         rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
         stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{P}{the parameters of the Narushin-Romanov-Griffin equation.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
}
\details{
The formula of the volume (\eqn{V}) of an egg based on the Narushin-Romanov-Griffin equation is:
\deqn{V(x)=\pi\int_{-A/2}^{A/2}y^2\,dx,}
where \eqn{y} denotes the Narushin-Romanov-Griffin equation (i.e., \code{\link{NRGE}}), and 
\eqn{A} denotes the egg's length, which is the first element in the parameter vector, \code{P}.  
}
\note{
The argument \code{P} in the \code{\link{VolumeNRGE}} function has the same parameters, as those in the 
\code{\link{NRGE}} function.  
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Narushin, V.G., Romanov, M.N., Griffin, D.K. (2021) Egg and math: introducing a universal formula for egg shape. 
\emph{Annals of the New York Academy of Sciences} 1505, 169\eqn{-}177. \doi{10.1111/nyas.14680}

Narushin, V.G., Romanov, M.N., Mishra, B., Griffin, D.K. (2022) Mathematical progression of 
avian egg shape with associated area and volume determinations. 
\emph{Annals of the New York Academy of Sciences} 1513, 65\eqn{-}78. \doi{10.1111/nyas.14771}

Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2023) 
A simple way to calculate the volume and surface area of avian eggs. 
\emph{Annals of the New York Academy of Sciences} 1524, 118\eqn{-}131. \doi{10.1111/nyas.15000}
}
\seealso{
\code{\link{curveNRGE}}, \code{\link{fitNRGE}}, \code{\link{NRGE}}, \code{\link{SurfaceAreaNRGE}}
}
\examples{
  Par6 <- c(4.51, 3.18, 0.1227, 2.2284)
  VolumeNRGE(P=Par6)

  myfun <- function(x){
    pi*NRGE(P=Par6, x=x)^2
  }
  integrate(myfun, -4.51/2, 4.51/2)$value
}
