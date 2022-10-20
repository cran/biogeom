\name{MLRFE}
\encoding{UTF-8}
\alias{MLRFE}
\title{
Modified Lobry-Rosso-Flandrois (LRF) Equation
}
\description{
\code{MLRFE} is used to calculate \eqn{y} values at given \eqn{x} values 
using the modified LRF equation or one of its simplified versions.
}
\usage{
MLRFE(P, x, simpver = 1)
}
\arguments{
  \item{P}{the parameters of the modified LRF equation or one of its simplified versions.}
  \item{x}{the given \eqn{x} values.}
  \item{simpver}{an optional argument to use the simplified version of the modified LRF equation.}
}
\details{

When \code{simpver = NULL}, the modified LRF equation is selected: 

\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, \frac{x_{\mathrm{min}}+x_{\mathrm{max}}}{2}\right)},}
\deqn{y = y_{\mathrm{opt}}\left\{\frac{\left(x-x_{\mathrm{min}}\right)\left(x-x_{\mathrm{max}}\right)^{2}}{\left(x_{\mathrm{max}}-x_{\mathrm{opt}}\right)\left[(x_{\mathrm{max}}-x_{\mathrm{opt}})(x-x_{\mathrm{opt}})-(x_{\mathrm{min}}-x_{\mathrm{opt}})(x_{\mathrm{opt}}+x_{\mathrm{max}}-2x)\right]}\right\}^{\delta};}
\deqn{\mbox{if } x \in{\left[\frac{x_{\mathrm{min}}+x_{\mathrm{max}}}{2}, x_{\mathrm{max}}\right)},}
\deqn{y = y_{\mathrm{opt}}\left\{\frac{\left(x-x_{\mathrm{max}}\right)\left(x-x_{\mathrm{min}}\right)^{2}}{\left(x_{\mathrm{opt}}-x_{\mathrm{min}}\right)\left[(x_{\mathrm{opt}}-x_{\mathrm{min}})(x-x_{\mathrm{opt}})-(x_{\mathrm{opt}}-x_{\mathrm{max}})(x_{\mathrm{opt}}+x_{\mathrm{min}}-2x)\right]}\right\}^{\delta};}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}

Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
\eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{min}}}, and \eqn{x_{\mathrm{max}}} are constants to be estimated;  
\eqn{y_{\mathrm{opt}}} represents the maximum \eqn{y}, and \eqn{x_{\mathrm{opt}}} is the \eqn{x} value associated with 
the maximum \eqn{y} (i.e., \eqn{y_{\mathrm{opt}}}); 
and \eqn{x_{\mathrm{min}}} and \eqn{x_{\mathrm{max}}} represents the 
lower and upper intersections between the curve and the \eqn{x}-axis. 
There are five elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{min}}}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta}, respectively.


When \code{simpver = 1}, the simplified version 1 is selected:

\deqn{\mbox{if } x \in{\left(0, \frac{x_{\mathrm{max}}}{2}\right)},}
\deqn{y = y_{\mathrm{opt}}\left\{\frac{x\left(x-x_{\mathrm{max}}\right)^{2}}{\left(x_{\mathrm{max}}-x_{\mathrm{opt}}\right)\left[(x_{\mathrm{max}}-x_{\mathrm{opt}})(x-x_{\mathrm{opt}})+x_{\mathrm{opt}}(x_{\mathrm{opt}}+x_{\mathrm{max}}-2x)\right]}\right\}^{\delta};}
\deqn{\mbox{if } x \in{\left[\frac{x_{\mathrm{max}}}{2}, x_{\mathrm{max}}\right)},}
\deqn{y = y_{\mathrm{opt}}\left\{\frac{\left(x-x_{\mathrm{max}}\right)x^{2}}{x_{\mathrm{opt}}\left[x_{\mathrm{opt}}(x-x_{\mathrm{opt}})-(x_{\mathrm{opt}}-x_{\mathrm{max}})(x_{\mathrm{opt}}-2x)\right]}\right\}^{\delta};}
\deqn{\mbox{if } x \notin{\left(0, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta}, respectively.


When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, \frac{x_{\mathrm{min}}+x_{\mathrm{max}}}{2}\right)},}
\deqn{y = \frac{y_{\mathrm{opt}}\left(x-x_{\mathrm{min}}\right)\left(x-x_{\mathrm{max}}\right)^{2}}{\left(x_{\mathrm{max}}-x_{\mathrm{opt}}\right)\left[(x_{\mathrm{max}}-x_{\mathrm{opt}})(x-x_{\mathrm{opt}})-(x_{\mathrm{min}}-x_{\mathrm{opt}})(x_{\mathrm{opt}}+x_{\mathrm{max}}-2x)\right]};}
\deqn{\mbox{if } x \in{\left[\frac{x_{\mathrm{min}}+x_{\mathrm{max}}}{2}, x_{\mathrm{max}}\right)},}
\deqn{y = \frac{y_{\mathrm{opt}}\left(x-x_{\mathrm{max}}\right)\left(x-x_{\mathrm{min}}\right)^{2}}{\left(x_{\mathrm{opt}}-x_{\mathrm{min}}\right)\left[(x_{\mathrm{opt}}-x_{\mathrm{min}})(x-x_{\mathrm{opt}})-(x_{\mathrm{opt}}-x_{\mathrm{max}})(x_{\mathrm{opt}}+x_{\mathrm{min}}-2x)\right]};}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, \eqn{x_{\mathrm{min}}}, and \eqn{x_{\mathrm{max}}}, respectively.


When \code{simpver = 3}, the simplified version 3 is selected:
\deqn{\mbox{if } x \in{\left(0, \frac{x_{\mathrm{max}}}{2}\right)},}
\deqn{y = \frac{y_{\mathrm{opt}}x\left(x-x_{\mathrm{max}}\right)^{2}}{\left(x_{\mathrm{max}}-x_{\mathrm{opt}}\right)\left[(x_{\mathrm{max}}-x_{\mathrm{opt}})(x-x_{\mathrm{opt}})+x_{\mathrm{opt}}(x_{\mathrm{opt}}+x_{\mathrm{max}}-2x)\right]};}
\deqn{\mbox{if } x \in{\left[\frac{x_{\mathrm{max}}}{2}, x_{\mathrm{max}}\right)},}
\deqn{y = \frac{y_{\mathrm{opt}}\left(x-x_{\mathrm{max}}\right)x^{2}}{x_{\mathrm{opt}}\left[x_{\mathrm{opt}}(x-x_{\mathrm{opt}})-(x_{\mathrm{opt}}-x_{\mathrm{max}})(x_{\mathrm{opt}}-2x)\right]};}
\deqn{\mbox{if } x \notin{\left(0, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are three elements in \code{P}, representing 
the values of \eqn{y_{\mathrm{opt}}}, \eqn{x_{\mathrm{opt}}}, and \eqn{x_{\mathrm{max}}}, respectively.
}
\value{
The \eqn{y} values predicted by the modified LRF equation or one of its simplified versions.
}
\note{
We have added n parameter \eqn{\delta} in the original LRF equation (i.e., \code{simpver = 2}) to increase the flexibility for data fitting.
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
\seealso{
\code{\link{areaovate}}, \code{\link{curveovate}}, \code{\link{fitovate}}, \code{\link{MbetaE}}, \code{\link{MBriereE}}, \code{\link{sigmoid}}
}
\examples{
x3   <- seq(-5, 15, len=2000)
Par3 <- c(3, 3, 10, 2)
y3   <- MbetaE(P=Par3, x=x3, simpver=1)

dev.new()
plot( x3, y3,cex.lab=1.5, cex.axis=1.5, type="l",
      xlab=expression(italic(x)), ylab=expression(italic(y)) )
 

}