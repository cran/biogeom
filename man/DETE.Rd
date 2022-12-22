\name{DETE}
\encoding{UTF-8}
\alias{DETE}
\title{
Calculation of the First-Order Derivative of the Explicit Troscianko Equation  
}
\description{
\code{DETE} is used to calculate the first-order derivative of the explicit Troscianko equation at a given \emph{x}-value.
}
\usage{
DETE(P, x)
}
\arguments{
  \item{P}{the parameters of the explicit Troscianko equation.}
  \item{x}{the \emph{x}-value used in the explicit Troscianko equation.}
}
\details{
The first-order derivative of the explicit Troscianko equation at a given \emph{x}-value is: 

\deqn{h(x) = \left\{\alpha_{1}+\frac{2\,\alpha_{2}}{a}\,x-\frac{x}{a}\,\left[1-\left(\frac{x}{a}\right)^2\right]^{-1}\right\}\,\exp\left\{\alpha_{0}+\alpha_{1}\left(\frac{x}{a}\right)+\alpha_{2}\left(\frac{x}{a}\right)^2\right\}\,\sqrt{ 1-\left(\frac{x}{a}\right)^2 }, }

where \code{P} has five parameters: \eqn{a}, \eqn{\alpha_{0}}, \eqn{\alpha_{1}}, and \eqn{\alpha_{2}}.
}
\note{
The argument \code{P} in the \code{\link{DETE}} function has the same parameters, as those in the 
\code{\link{ETE}} function.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2022) 
Can we simply calculate the volume and surface area of an egg? Evidence from six species of poultry. 
\emph{Annals of the New York Academy of Sciences} Under preparation.
}
\seealso{
\code{\link{ETE}}, \code{\link{fitETE}}, \code{\link{SurfaceAreaETE}}
}
\examples{
  Par5 <- c(2.25, -0.38, -0.29, -0.16)
  xx2  <- seq(-2.25, 2.25, by=0.001)
  h1   <- DETE(P=Par5, x=xx2)
  h2   <- -DETE(P=Par5, x=xx2)
  ind  <- which(is.na(h1) | is.na(h2))
  xx2  <- xx2[-ind]
  h1   <- h1[-ind]
  h2  <- h2[-ind]

  dev.new()
  plot(xx2, h1, type="l", col=4, cex.lab=1.5, cex.axis=1.5,
       xlim=c(-2.25, 2.25), ylim=c(-30, 30), xlab=expression(italic(x)), 
       ylab=expression(paste(italic(h), "(", italic(x), ")", sep="")))
  lines(xx2, h2, col=2)  

  graphics.off()
}