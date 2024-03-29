\name{SCSE}
\encoding{UTF-8}
\alias{SCSE}
\title{
Sarabia-Castillo-Slottje Equation (SCSE)
}
\description{
\code{SCSE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the Sarabia-Castillo-Slottje equation. The equation describes the \eqn{y} coordinates of the Lorenz curve.
}
\usage{
SCSE(P, x)
}
\arguments{
  \item{P}{the parameters of the Sarabia-Castillo-Slottje equation.}
  \item{x}{the given \eqn{x} values ranging between 0 and 1.}
}
\details{

\deqn{y = x^{\gamma}\left[1-\left(1-x\right)^{\alpha}\right]^{\beta}.}

Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
and \eqn{\gamma}, \eqn{\alpha} and \eqn{\beta} are constants to be estimated, where \eqn{\gamma \ge 0}, 
\eqn{0 < \alpha \le 1}, and \eqn{\beta \ge 1}.   
There are three elements in \code{P}, representing 
the values of \eqn{\gamma}, \eqn{\alpha} and \eqn{\beta}, respectively.
}
\value{
The \eqn{y} values predicted by the Sarabia-Castillo-Slottje equation.
}
\note{
The numerical range of \eqn{x} should range between 0 and 1.
}
\references{
Sarabia, J.-M., Castillo, E., Slottje, D.J. (1999) An ordered family of Lorenz curves. 
\emph{Journal of Econometrics}. 91, 43\eqn{-}60. \doi{10.1016/S0304-4076(98)00048-7}

Sitthiyot, T., Holasut, K. (2023) A universal model for the Lorenz curve with novel applications 
for datasets containing zeros and/or exhibiting extreme inequality. \emph{Scientific Reports} 
13, 4729. \doi{10.1038/s41598-023-31827-x}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\seealso{
\code{\link{fitLorenz}}, \code{\link{MPerformanceE}}, \code{\link{SarabiaE}}, \code{\link{SHE}}
}
\examples{
X1  <- seq(0, 1, len=2000)
Pa2 <- c(0, 0.790, 1.343)
Y2  <- SCSE(P=Pa2, x=X1)

dev.new()
plot( X1, Y2, cex.lab=1.5, cex.axis=1.5, type="l", asp=1, xaxs="i", 
      yaxs="i", xlim=c(0, 1), ylim=c(0, 1), 
      xlab="Cumulative proportion of the number of infructescences", 
      ylab="Cumulative proportion of the infructescence length" )

graphics.off()
}