\name{MbetaE}
\encoding{UTF-8}
\alias{MbetaE}
\title{
Modified Beta Equation
}
\description{
\code{MbetaE} is used to calculate \eqn{y} values at given \eqn{x} values 
using the modified beta equation or one of its simplified versions.
}
\usage{
MbetaE(P, x, simpver = 1)
}
\arguments{
  \item{P}{the parameters of the modified beta equation or one of its simplified versions.}
  \item{x}{the given \eqn{x} values.}
  \item{simpver}{an optional argument to use the simplified version of the modified beta equation.}
}
\details{

When \code{simpver = NULL}, the modified beta equation is selected: 
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = y_{\mathrm{opt}}{ \left[\left(\frac{x_{\mathrm{max}}-x}{x_{\mathrm{max}}-x_{\mathrm{opt}}}\right)\left(\frac{x-x_{\mathrm{min}}}{x_{\mathrm{opt}}-x_{\mathrm{min}}}\right)^{\frac{x_{\mathrm{opt}}-x_{\mathrm{min}}}{x_{\mathrm{max}}-x_{\mathrm{opt}}}} \right] }^{\delta};}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
\eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{min}}}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta} are constants to be estimated;  
\eqn{y_{\mathrm{opt}}} represents the maximum \eqn{y}, and \eqn{x_{\mathrm{opt}}} is the \eqn{x} value associated with 
the maximum \eqn{y} (i.e., \eqn{y_{\mathrm{opt}}}); 
and \eqn{x_{\mathrm{min}}} and \eqn{x_{\mathrm{max}}} represent the 
lower and upper intersections between the curve and the \eqn{x}-axis. \eqn{y} is defined as 0 
when \eqn{x < x_{\mathrm{min}}} or \eqn{x > x_{\mathrm{max}}}. There are five elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{min}}}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta}, respectively.

\eqn{\quad} When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{\mbox{if } x \in{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = y_{\mathrm{opt}}{ \left[\left(\frac{x_{\mathrm{max}}-x}{x_{\mathrm{max}}-x_{\mathrm{opt}}}\right)\left(\frac{x}{x_{\mathrm{opt}}}\right)^{\frac{x_{\mathrm{opt}}}{x_{\mathrm{max}}-x_{\mathrm{opt}}}} \right] }^{\delta};}
\deqn{\mbox{if } x \notin{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta}, respectively.


\eqn{\quad} When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = y_{\mathrm{opt}}{ \left(\frac{x_{\mathrm{max}}-x}{x_{\mathrm{max}}-x_{\mathrm{opt}}}\right)\left(\frac{x-x_{\mathrm{min}}}{x_{\mathrm{opt}}-x_{\mathrm{min}}}\right)^{\frac{x_{\mathrm{opt}}-x_{\mathrm{min}}}{x_{\mathrm{max}}-x_{\mathrm{opt}}}} };}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{min}}}, and \eqn{x_{\mathrm{max}}}, respectively.


\eqn{\quad} When \code{simpver = 3}, the simplified version 3 is selected:
\deqn{\mbox{if } x \in{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = y_{\mathrm{opt}}{ \left(\frac{x_{\mathrm{max}}-x}{x_{\mathrm{max}}-x_{\mathrm{opt}}}\right)\left(\frac{x}{x_{\mathrm{opt}}}\right)^{\frac{x_{\mathrm{opt}}}{x_{\mathrm{max}}-x_{\mathrm{opt}}}} }; } 
\deqn{\mbox{if } x \notin{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are three elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, and \eqn{x_{\mathrm{max}}}, respectively.
}
\value{
The \eqn{y} values predicted by the modified beta equation or one of its simplified versions.
}
\references{
Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\note{
We have added a parameter \eqn{\delta} in the original beta equation (i.e., \code{simpver = 2}) to increase the flexibility for data fitting.
}
\seealso{
\code{\link{areaovate}}, \code{\link{curveovate}}, \code{\link{fitovate}}, \code{\link{fitsigmoid}}, 
\code{\link{MBriereE}}, \code{\link{MLRFE}}, \code{\link{MPerformanceE}}, \code{\link{sigmoid}}
}
\examples{
x1   <- seq(-5, 15, len=2000)
Par1 <- c(3, 3, 10, 2)
y1   <- MbetaE(P=Par1, x=x1, simpver=1)

dev.new()
plot( x1, y1,cex.lab=1.5, cex.axis=1.5, type="l",
      xlab=expression(italic(x)), ylab=expression(italic(y)) )
 
graphics.off()
}