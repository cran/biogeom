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
  \item{P}{the parameters of the Preston equation or one of its simplified versions.}
  \item{x}{the \emph{x}-value used in the explicit Preston equation.}
  \item{simpver}{an optional argument to use the simplified version of the explicit Preston equation.}
}
\details{
When \code{simpver = NULL}, the first-order derivative of 
the explicit Preston equation at a given \emph{x}-value is selected: 
\deqn{ f(x)=\frac{b\left[a^{4}\,c_{1}+a^{3}\left(2\,c_{2}-1\right)x+
a^2\left(3\,c_{3}-2\,c_{1}\right)x^{2}-3\,a\,c_{2}x^3-4\,c_{3}\,x^{4}\right]}{a^4\sqrt{a^2-x^2}} }
where \code{P} has five parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, \eqn{c_{2}}, and \eqn{c_{3}}.

\eqn{\quad} When \code{simpver = 1}, the first-order derivative of the simplified version 1 is selected: 
\deqn{ f(x)=\frac{b\left[a^{4}\,c_{1}+a^{3}\left(2\,c_{2}-1\right)x-
2\,a^2\,c_{1}\,x^{2}-3\,a\,c_{2}x^3\right]}{a^4\sqrt{a^2-x^2}} }
where \code{P} has four parameters: \eqn{a}, \eqn{b}, \eqn{c_{1}}, and \eqn{c_{2}}.

\eqn{\quad} When \code{simpver = 2}, the first-order derivative of the simplified version 2 is selected: 
\deqn{ f(x)=\frac{b\left[a^{4}\,c_{1}-a^{3}\,x-
2\,a^2\,c_{1}\,x^{2}\right]}{a^4\sqrt{a^2-x^2}} }
where \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{1}}.

\eqn{\quad} When \code{simpver = 3}, the first-order derivative of the simplified version 3 is selected: 
\deqn{ f(x)=\frac{b\left[a^{3}\left(2\,c_{2}-1\right)x-3\,a\,c_{2}x^3\right]}{a^4\sqrt{a^2-x^2}} }
where \code{P} has three parameters: \eqn{a}, \eqn{b}, and \eqn{c_{2}}.
}
\note{
The \code{\link{DEPE}} function has the same parameters, \code{P}, as those in the 
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
\code{\link{EPE}}, \code{\link{PE}}, \code{\link{SurfaceAreaPE}}, \code{\link{VolumePE}}
}
\examples{
  Par3 <- c(4.27, 2.90, 0.0868, 0.0224, -0.0287)
  x    <- seq(-4.27, 4.27, by=0.001)
  y3   <- DEPE(P=Par3, x=x, simpver=NULL)
  y4   <- -DEPE(P=Par3, x=x, simpver=NULL)

  dev.new()
  plot(x, y3, type="l", col=4, cex.lab=1.5, cex.axis=1.5,
       xlim=c(-5, 5), ylim=c(-35, 35), 
       xlab=expression(italic(x)), ylab=expression(italic(f)))
  lines(x, y4, col=2)  
}