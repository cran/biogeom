\name{MBriereE}
\encoding{UTF-8}
\alias{MBriereE}
\title{
Modified Briere Equation
}
\description{
\code{MBriereE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the modified Brière equation or one of its simplified versions. 
}
\usage{
MBriereE(P, x, simpver = 1)
}
\arguments{
  \item{P}{the parameters of the modified Brière equation or one of its simplified versions.}
  \item{x}{the given \eqn{x} values.}
  \item{simpver}{an optional argument to use the simplified version of the modified Brière equation.}
}
\details{

When \code{simpver = NULL}, the modified Brière equation is selected: 
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, x_{\mathrm{max}}\right)},}
\deqn{y = a\left|x(x-x_{\mathrm{min}})(x_{\mathrm{max}}-x)^{1/m}\right|^{\delta};}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
and \eqn{a}, \eqn{m}, \eqn{x_{\mathrm{min}}}, and \eqn{x_{\mathrm{max}}} are constants to be estimated, 
where \eqn{x_{\mathrm{min}}} and \eqn{x_{\mathrm{max}}} represents the 
lower and upper intersections between the curve and the \eqn{x}-axis. \eqn{y} is defined as 0 
when \eqn{x < x_{\mathrm{min}}} or \eqn{x > x_{\mathrm{max}}}. There are five elements in \code{P}, representing 
the values of \eqn{a}, \eqn{m}, \eqn{x_{\mathrm{min}}}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta}, respectively.

When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{\mbox{if } x \in{\left(0, x_{\mathrm{max}}\right)},}
\deqn{y = a\left|x^{2}(x_{\mathrm{max}}-x)^{1/m}\right|^{\delta};}
\deqn{\mbox{if } x \notin{\left(0, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P}, representing 
the values of \eqn{a}, \eqn{m}, \eqn{x_{\mathrm{max}}}, and \eqn{\delta}, respectively.


When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, x_{\mathrm{max}}\right)},}
\deqn{y = ax(x-x_{\mathrm{min}})(x_{\mathrm{max}}-x)^{1/m};}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P} representing 
the values of \eqn{a}, \eqn{m}, \eqn{x_{\mathrm{min}}}, and \eqn{x_{\mathrm{max}}}, respectively.


When \code{simpver = 3}, the simplified version 3 is selected: 
\deqn{\mbox{if } x \in{\left(0, x_{\mathrm{max}}\right)},}
\deqn{y = ax^{2}(x_{\mathrm{max}}-x)^{1/m};}
\deqn{\mbox{if } x \notin{\left(0, x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are three elements in \code{P} representing 
the values of \eqn{a}, \eqn{m}, and \eqn{x_{\mathrm{max}}}, respectively.
}
\value{
The \eqn{y} values predicted by the modified Brière equation or one of its simplified versions.
}
\note{
We have added a parameter \eqn{\delta} in the original Brière equation (i.e., \code{simpver = 2}) to increase the flexibility for data fitting.
}
\references{
Brière, J.-F., Pracros, P, Le Roux, A.-Y., Pierre, J.-S. (1999) A novel rate
model of temperature-dependent development for arthropods. \emph{Environmental 
Entomology} 28, 22\eqn{-}29. \doi{10.1093/ee/28.1.22}

Cao, L., Shi, P., Li, L., Chen, G. (2019) A new flexible sigmoidal growth model. \emph{Symmetry} 11, 204.
\doi{10.3390/sym11020204}

Jin, J., Quinn, B.K., Shi, P. (2022) The modified Brière equation and its 
applications. \emph{Plants} 11, 1769. \doi{10.3390/plants11131769}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\seealso{
\code{\link{areaovate}}, \code{\link{curveovate}}, \code{\link{fitovate}}, \code{\link{MbetaE}}, \code{\link{MLRFE}}, \code{\link{sigmoid}}
}
\examples{
x2   <- seq(-5, 15, len=2000)
Par2 <- c(0.01, 3, 0, 10, 1)
y2   <- MBriereE(P=Par2, x=x2, simpver=NULL)

dev.new()
plot( x2, y2, cex.lab=1.5, cex.axis=1.5, type="l",
      xlab=expression(italic(x)), ylab=expression(italic(y)) )

graphics.off()
}