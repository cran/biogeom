\name{SHE}
\encoding{UTF-8}
\alias{SHE}
\title{
Sitthiyot-Holasut Equation
}
\description{
\code{SHE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the Sitthiyot-Holasut equation. The equation describes the \eqn{y} coordinates of the Lorenz curve.
}
\usage{
SHE(P, x)
}
\arguments{
  \item{P}{the parameters of the Sitthiyot-Holasut equation.}
  \item{x}{the given \eqn{x} values ranging between 0 and 1.}
}
\details{

\deqn{\mbox{if } x > \delta,}
\deqn{y = \left(1-\rho\right)\,\left[\left(\frac{2}{P+1}\right)\left(\frac{x-\delta}{1-\delta}\right)\right] + 
\rho\,\left[\left(1-\omega\right)\left(\frac{x-\delta}{1-\delta}\right)^{P}+\omega\,\left\{1-\left[1-\left(\frac{x-\delta}{1-\delta}\right)\right]^{\frac{1}{P}}\right\}\right];}
\deqn{\mbox{if } x \le \delta,}
\deqn{y = 0.}
Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
and \eqn{\delta}, \eqn{\rho}, \eqn{\omega} and \eqn{P} are constants to be estimated, where \eqn{0 \le \delta < 1}, 
\eqn{0 \le \rho \le 1}, \eqn{0 \le \omega \le 1}, and \eqn{P \ge 1}.   
There are four elements in \code{P}, representing 
the values of \eqn{\delta}, \eqn{\rho}, \eqn{\omega} and \eqn{P}, respectively.

}
\value{
The \eqn{y} values predicted by the Sitthiyot-Holasut equation.
}
\note{
The numerical range of \eqn{x} should range between 0 and 1. 
When \eqn{x < \delta}, the \eqn{x} value is assigned to be \eqn{\delta}. 
}
\references{
Sitthiyot, T., Holasut, K. (2023) A universal model for the Lorenz curve with novel applications 
for datasets containing zeros and/or exhibiting extreme inequality. \emph{Scientific Reports} 
13, 4729. \doi{10.1038/s41598-023-31827-x}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\seealso{
\code{\link{fitLorenz}}, \code{\link{MPerformanceE}}, \code{\link{SarabiaE}}, \code{\link{SCSE}}
}
\examples{
X1  <- seq(0, 1, len=2000)
Pa3 <- c(0, 1, 0.446, 1.739)
Y3  <- SHE(P=Pa3, x=X1)

dev.new()
plot( X1, Y3, cex.lab=1.5, cex.axis=1.5, type="l", asp=1, xaxs="i", 
      yaxs="i", xlim=c(0, 1), ylim=c(0, 1), 
      xlab="Cumulative proportion of the number of infructescences", 
      ylab="Cumulative proportion of the infructescence length" )

graphics.off()
}