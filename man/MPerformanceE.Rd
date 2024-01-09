\name{MPerformanceE}
\encoding{UTF-8}
\alias{MPerformanceE}
\title{
Modified Performance Equation
}
\description{
\code{MPerformanceE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the modified performance equation or one of its simplified versions. 
}
\usage{
MPerformanceE(P, x, simpver = 1)
}
\arguments{
  \item{P}{the parameters of the modified performance equation or one of its simplified versions.}
  \item{x}{the given \eqn{x} values.}
  \item{simpver}{an optional argument to use the simplified version of the modified performance equation.}
}
\details{

When \code{simpver = NULL}, the modified performance equation is selected: 
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = c\left(1-e^{-K_{1}\left(x-x_{\mathrm{min}}\right)}\right)^{a}\left(1-e^{K_{2}\left(x-x_{\mathrm{max}}\right)}\right)^{b};}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
Here, \eqn{x} and \eqn{y} represent the independent and dependent variables, respectively; 
and \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, \eqn{x_{\mathrm{min}}}, \eqn{x_{\mathrm{max}}}, \eqn{a}, 
and \eqn{b} are constants to be estimated, 
where \eqn{x_{\mathrm{min}}} and \eqn{x_{\mathrm{max}}} represents the 
lower and upper intersections between the curve and the \eqn{x}-axis. \eqn{y} is defined as 0 
when \eqn{x < x_{\mathrm{min}}} or \eqn{x > x_{\mathrm{max}}}. There are seven elements in \code{P}, representing 
the values of \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, \eqn{x_{\mathrm{min}}}, \eqn{x_{\mathrm{max}}}, \eqn{a}, and \eqn{b}, respectively.

\eqn{\quad} When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{\mbox{if } x \in{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = c\left(1-e^{-K_{1}x}\right)^{a}\left(1-e^{K_{2}\left(x-x_{\mathrm{max}}\right)}\right)^{b};}
\deqn{\mbox{if } x \notin{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are six elements in \code{P}, representing 
the values of \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, \eqn{x_{\mathrm{max}}}, \eqn{a}, and \eqn{b} respectively.


\eqn{\quad} When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{\mbox{if } x \in{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = c\left(1-e^{-K_{1}\left(x-x_{\mathrm{min}}\right)}\right)\left(1-e^{K_{2}\left(x-x_{\mathrm{max}}\right)}\right);}
\deqn{\mbox{if } x \notin{\left(x_{\mathrm{min}}, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are five elements in \code{P} representing 
the values of \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, \eqn{x_{\mathrm{min}}}, and \eqn{x_{\mathrm{max}}}, respectively.


\eqn{\quad} When \code{simpver = 3}, the simplified version 3 is selected: 
\deqn{\mbox{if } x \in{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = c\left(1-e^{-K_{1}x}\right)\left(1-e^{K_{2}\left(x-x_{\mathrm{max}}\right)}\right);}
\deqn{\mbox{if } x \notin{\left(0, \ x_{\mathrm{max}}\right)},}
\deqn{y = 0.}
There are four elements in \code{P} representing 
the values of \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, and \eqn{x_{\mathrm{max}}}, respectively.

\eqn{\quad} When \code{simpver = 4}, the simplified version 4 is selected: 
\deqn{\mbox{if } x \in{\left(0, \ \sqrt{2}\right)},}
\deqn{y = c\left(1-e^{-K_{1}x}\right)^{a}\left(1-e^{K_{2}\left(x-\sqrt{2}\right)}\right)^{b};}
\deqn{\mbox{if } x \notin{\left(0, \ \sqrt{2}\right)},}
\deqn{y = 0.}
There are five elements in \code{P}, representing 
the values of \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, \eqn{a}, and \eqn{b}, respectively.

\eqn{\quad} When \code{simpver = 5}, the simplified version 5 is selected:
\deqn{\mbox{if } x \in{\left(0, \ \sqrt{2}\right)},}
\deqn{y = c\left(1-e^{-K_{1}x}\right)\left(1-e^{K_{2}\left(x-\sqrt{2}\right)}\right);}
\deqn{\mbox{if } x \notin{\left(0, \ \sqrt{2}\right)},}
\deqn{y = 0.}
There are three elements in \code{P}, representing 
the values of \eqn{c}, \eqn{K_{1}}, and \eqn{K_{2}}, respectively.

}
\value{
The \eqn{y} values predicted by the modified performance equation or one of its simplified versions.
}
\note{
We have added two parameters \eqn{a} and \eqn{b} in the original performance 
equation (i.e., \code{simpver = 2}) to increase the flexibility for data fitting. 
The cases of \code{simpver = 4} and \code{simpver = 5} are used to describe the rotated and right-shifted 
Lorenz curve (see Lian et al. [2023] for details).
}
\references{
Huey, R.B., Stevenson, R.D. (1979) Integrating thermal physiology and ecology of ectotherms: 
a discussion of approaches. \emph{American Zoologist} 19, 357\eqn{-}366. \doi{10.1093/icb/19.1.357}

Lian, M., Shi, P., Zhang, L., Yao, W., Gielis, J., Niklas, K.J. (2023) A generalized performance equation 
and its application in measuring the Gini index of leaf size inequality. 
\emph{Trees \eqn{-} Structure and Function} 37, 1555\eqn{-}1565. \doi{10.1007/s00468-023-02448-8}

Shi, P., Ge, F., Sun, Y., Chen, C. (2011) A simple model for describing 
the effect of temperature on insect developmental rate. \emph{Journal of Asia-Pacific Entomology} 
14, 15\eqn{-}20. \doi{10.1016/j.aspen.2010.11.008} 

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\seealso{
\code{\link{areaovate}}, \code{\link{curveovate}}, \code{\link{fitLorenz}}, 
\code{\link{fitovate}}, \code{\link{fitsigmoid}}, \code{\link{MbetaE}}, 
\code{\link{MBriereE}}, \code{\link{MLRFE}}, \code{\link{sigmoid}}
}
\examples{
x4   <- seq(0, 40, len=2000)
Par4 <- c(0.117, 0.090, 0.255, 5, 35, 1, 1)
y4   <- MPerformanceE(P=Par4, x=x4, simpver=NULL)

dev.new()
plot( x4, y4, cex.lab=1.5, cex.axis=1.5, type="l",
      xlab=expression(italic(x)), ylab=expression(italic(y)) )

graphics.off()
}