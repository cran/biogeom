\name{VolumeEPE}
\encoding{UTF-8}
\alias{VolumeEPE}
\title{
Calculation of the Volume of An Egg Based on the Explicit Preston Equation
}
\description{
\code{VolumeEPE} is used to calculate the volume of an egg that follows the explicit Preston equation.
}
\usage{
VolumeEPE(P, simpver = NULL)
}
\arguments{
  \item{P}{the parameters of the explicit Preston equation or one of its simplified versions.}
  \item{simpver}{an optional argument to use the simplified version of the explicit Preston equation.}
}
\details{
When \code{simpver = NULL}, the volume formula (\eqn{V}) of the explicit Preston equation is selected: 
\deqn{ V(x) = \frac{4\,\pi}{315}a\,b^{2}\left(105+21\,c_{1}^{2}+42\,c_{2}+9\,c_{2}^2+18\,c_{1}\,c_{3}+5\,c_{3}^2\right), }
where \code{P} has five parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, \eqn{c_{2}}, and \eqn{c_{3}}.

\eqn{\quad} When \code{simpver = 1}, the volume formula of the simplified version 1 is selected:
\deqn{ V(x) = \frac{4\,\pi}{315}a\,b^{2}\left(105+21\,c_{1}^{2}+42\,c_{2}+9\,c_{2}^2\right), }
where \code{P} has four parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, and \eqn{c_{2}}.

\eqn{\quad} When \code{simpver = 2}, the volume formula of the simplified version 2 is selected:
\deqn{ V(x) = \frac{4\,\pi}{315}a\,b^{2}\left(105+21\,c_{1}^{2}\right), }
where \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{1}}.

\eqn{\quad} When \code{simpver = 3}, the volume formula of the simplified version 3 is selected:
\deqn{ V(x) = \frac{4\,\pi}{315}a\,b^{2}\left(105+42\,c_{2}+9\,c_{2}^2\right), }
where \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{2}}.
}
\note{
The argument \code{P} in the \code{\link{VolumeEPE}} function has the same parameters, as those in the 
\code{\link{EPE}} function. 
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
\code{\link{EPE}}, \code{\link{fitEPE}}, \code{\link{SurfaceAreaEPE}}
}
\examples{
  Par3 <- c(4.27, 2.90, 0.0868, 0.0224, -0.0287)
  VolumeEPE(P=Par3, simpver=NULL)

  # Test the case when simpver = NULL
  a    <- Par3[1]
  b    <- Par3[2]
  c1   <- Par3[3]
  c2   <- Par3[4]
  c3   <- Par3[5]
  pi*4/315*a*b^2*(105+21*c1^2+42*c2+9*c2^2+18*c1*c3+5*c3^2)

  myfun <- function(x){
    pi*EPE(P=Par3, x=x, simpver=NULL)^2
  }
  integrate(myfun, -4.27, 4.27)$value
}
