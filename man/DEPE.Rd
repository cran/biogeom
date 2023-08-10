\name{DEPE}
\encoding{UTF-8}
\alias{DEPE}
\title{
Calculation of the First-Order Derivative of the Explicit Preston Equation  
}
\description{
\code{DEPE} is used to calculate the first-order derivative of the explicit Preston equation at a given \emph{x}-value.
}
\usage{
DEPE(P, x, simpver = NULL)
}
\arguments{
  \item{P}{the parameters of the explicit Preston equation or one of its simplified versions.}
  \item{x}{the \emph{x}-value used in the explicit Preston equation.}
  \item{simpver}{an optional argument to use the simplified version of the explicit Preston equation.}
}
\details{
When \code{simpver = NULL}, the first-order derivative of 
the explicit Preston equation at a given \emph{x}-value is selected: 
\deqn{ f(x)=\frac{b\left[a^{4}\,c_{1}+a^{3}\left(2\,c_{2}-1\right)x+
a^2\left(3\,c_{3}-2\,c_{1}\right)x^{2}-3\,a\,c_{2}x^3-4\,c_{3}\,x^{4}\right]}{a^4\sqrt{a^2-x^2}}, }
where \code{P} has five parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, \eqn{c_{2}}, and \eqn{c_{3}}.

\eqn{\quad} When \code{simpver = 1}, the first-order derivative of the simplified version 1 is selected: 
\deqn{ f(x)=\frac{b\left[a^{4}\,c_{1}+a^{3}\left(2\,c_{2}-1\right)x-
2\,a^2\,c_{1}\,x^{2}-3\,a\,c_{2}x^3\right]}{a^4\sqrt{a^2-x^2}}, }
where \code{P} has four parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, and \eqn{c_{2}}.

\eqn{\quad} When \code{simpver = 2}, the first-order derivative of the simplified version 2 is selected: 
\deqn{ f(x)=\frac{b\left[a^{4}\,c_{1}-a^{3}\,x-
2\,a^2\,c_{1}\,x^{2}\right]}{a^4\sqrt{a^2-x^2}}, }
where \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{1}}.

\eqn{\quad} When \code{simpver = 3}, the first-order derivative of the simplified version 3 is selected: 
\deqn{ f(x)=\frac{b\left[a^{3}\left(2\,c_{2}-1\right)x-3\,a\,c_{2}x^3\right]}{a^4\sqrt{a^2-x^2}}, }
where \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{2}}.
}
\note{
The argument \code{P} in the \code{\link{DEPE}} function has the same parameters, as those in the 
\code{\link{EPE}} function. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2023) 
A simple way to calculate the volume and surface area of avian eggs. 
\emph{Annals of the New York Academy of Sciences} 1524, 118\eqn{-}131. \doi{10.1111/nyas.15000}
}
\seealso{
\code{\link{EPE}}, \code{\link{fitEPE}}, \code{\link{SurfaceAreaEPE}}
}
\examples{
  Par3 <- c(4.27, 2.90, 0.0868, 0.0224, -0.0287)
  xx1  <- seq(-4.27, 4.27, by=0.001)
  f1   <- DEPE(P=Par3, x=xx1, simpver=NULL)
  f2   <- -DEPE(P=Par3, x=xx1, simpver=NULL)

  dev.new()
  plot(xx1, f1, type="l", col=4, cex.lab=1.5, cex.axis=1.5,
       xlim=c(-5, 5), ylim=c(-35, 35), xlab=expression(italic(x)), 
       ylab=expression(paste(italic(f), "(", italic(x), ")", sep="")))
  lines(xx1, f2, col=2)  

  graphics.off()
}