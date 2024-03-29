\name{SarabiaE}
\encoding{UTF-8}
\alias{SarabiaE}
\title{
Sarabia Equation
}
\description{
\code{SarabiaE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the Sarabia equation. The equation describes the \eqn{y} coordinates of the Lorenz curve.
}
\usage{
SarabiaE(P, x)
}
\arguments{
  \item{P}{the parameters of the Sarabia equation.}
  \item{x}{the given \eqn{x} values ranging between 0 and 1.}
}
\details{

\deqn{y = \left(1-\lambda+\eta\right)x+\lambda x^{a_1 + 1}-\eta \left[1-\left(1-x\right)^{a_2 + 1}\right].}

Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
and \eqn{\lambda}, \eqn{\eta}, \eqn{a_1} and \eqn{a_2} are constants to be estimated, where 
\eqn{a_1 \ge 0}, \eqn{a_2 + 1 \ge 0}, \eqn{\eta\,a_2 + \lambda \le 1}, \eqn{\lambda \ge 0}, and \eqn{\eta\,a_2 \ge 0}.   
There are four elements in \code{P}, representing 
the values of \eqn{\lambda}, \eqn{\eta}, \eqn{a_1} and \eqn{a_2}, respectively.

}
\value{
The \eqn{y} values predicted by the Sarabia equation.
}
\note{
The numerical range of \eqn{x} should range between 0 and 1. 
}
\references{
Sarabia, J.-M. (1997) A hierarchy of Lorenz curves based on the generalized Tukey's lambda distribution. 
\emph{Econometric Reviews} 16, 305\eqn{-}320. \doi{10.1080/07474939708800389} 

Sitthiyot, T., Holasut, K. (2023) A universal model for the Lorenz curve with novel applications 
for datasets containing zeros and/or exhibiting extreme inequality. \emph{Scientific Reports} 
13, 4729. \doi{10.1038/s41598-023-31827-x}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\seealso{
\code{\link{fitLorenz}}, \code{\link{MPerformanceE}}, \code{\link{SCSE}}, \code{\link{SHE}}
}
\examples{
X1  <- seq(0, 1, len=2000)
Pa1 <- c(0.295, 101.485, 0.705, 0.003762)
Y1  <- SarabiaE(P=Pa1, x=X1)

dev.new()
plot( X1, Y1, cex.lab=1.5, cex.axis=1.5, type="l", asp=1, xaxs="i", 
      yaxs="i", xlim=c(0, 1), ylim=c(0, 1), 
      xlab="Cumulative proportion of the number of infructescences", 
      ylab="Cumulative proportion of the infructescence length" )

graphics.off()
}