\name{PE}
\encoding{UTF-8}
\alias{PE}
\title{
Calculation of the Abscissa, Ordinate and Distance From the Origin For an Arbitrary Point on the Preston Curve 
}
\description{
\code{PE} is used to calculate the abscissa, ordinate and distance from the origin for an arbitrary point on 
the Preston curve that was generated by the original Preston equation or one of its simplified versions at a given angle.
}
\usage{
PE(P, zeta, simpver = NULL)
}
\arguments{
  \item{P}{the parameters of the original Preston equation or one of its simplified versions.}
  \item{zeta}{the angle(s) used in the Preston equation.}
  \item{simpver}{an optional argument to use the simplified version of the original Preston equation.}
}
\details{
When \code{simpver = NULL}, the original Preston equation is selected: 
\deqn{y = a\ \mathrm{sin}\,\zeta, }
\deqn{x = b\ \mathrm{cos}\,\zeta\left(1+c_{1}\,\mathrm{sin}\,\zeta+c_{2}\,\mathrm{sin}^{2}\,\zeta+c_{3}\,\mathrm{sin}^{3}\,\zeta\right), }
\deqn{r = \sqrt{x^{2}+y^{2}}, }
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Preston curve  
corresponding to an angle \eqn{\zeta}; \eqn{r} represents the distance of the point from the origin; \eqn{a}, \eqn{b}, \eqn{c_{1}}, 
\eqn{c_{2}}, and \eqn{c_{3}} are parameters to be estimated.

\eqn{\quad} When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{y = a\ \mathrm{sin}\,\zeta, }
\deqn{x = b\ \mathrm{cos}\,\zeta\left(1+c_{1}\,\mathrm{sin}\,\zeta+c_{2}\,\mathrm{sin}^{2}\,\zeta\right), }
\deqn{r = \sqrt{x^{2}+y^{2}}, }
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Preston curve  
corresponding to an angle \eqn{\zeta}; \eqn{r} represents the distance of the point from the origin; \eqn{a}, \eqn{b}, \eqn{c_{1}}, 
and \eqn{c_{2}} are parameters to be estimated.

\eqn{\quad} When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{y = a\ \mathrm{sin}\,\zeta, }
\deqn{x = b\ \mathrm{cos}\,\zeta\left(1+c_{1}\,\mathrm{sin}\,\zeta\right), }
\deqn{r = \sqrt{x^{2}+y^{2}}, }
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Preston curve  
corresponding to an angle \eqn{\zeta}; \eqn{r} represents the distance of the point from the origin; \eqn{a}, \eqn{b}, and \eqn{c_{1}}  
are parameters to be estimated.

\eqn{\quad} When \code{simpver = 3}, the simplified version 3 is selected:
\deqn{y = a\ \mathrm{sin}\,\zeta, }
\deqn{x = b\ \mathrm{cos}\,\zeta\left(1+c_{2}\,\mathrm{sin}^{2}\,\zeta\right), }
\deqn{r = \sqrt{x^{2}+y^{2}}, }
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Preston curve  
corresponding to an angle \eqn{\zeta}; \eqn{r} represents the distance of the point from the origin; \eqn{a}, \eqn{b}, and 
\eqn{c_{2}} are parameters to be estimated.

}
\value{
  \item{x}{the abscissa(s) of the Preston curve corresponding to the given angle(s).}
  \item{y}{the ordinate(s) of the Preston curve corresponding to the given angle(s).}
  \item{r}{the distance(s) of the Preston curve corresponding to the given angle(s) from the origin.}
}
\note{
\eqn{\zeta} is NOT the polar angle corresponding to \eqn{r}, i.e., 
\deqn{y \neq r\,\mathrm{sin}\,\zeta, } 
\deqn{x \neq r\,\mathrm{cos}\,\zeta. } 
Let \eqn{\varphi} be the polar angle corresponding to \eqn{r}. We have: 
\deqn{ \zeta = \mathrm{arc\,sin}\frac{ r\ \mathrm{sin}\,\varphi }{a}. }
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Biggins, J.D., Montgomeries, R.M., Thompson, J.E., Birkhead, T.R. (2022) 
Preston's universal formula for avian egg shape. \emph{Ornithology} 
139, ukac028. \doi{10.1093/ornithology/ukac028}

Biggins, J.D., Thompson, J.E., Birkhead, T.R. (2018) Accurately quantifying 
the shape of birds' eggs. \emph{Ecology and Evolution} 8, 9728\eqn{-}9738. \doi{10.1002/ece3.4412}

Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Wang, L., Quinn, B.K., Gielis, J. (2023) A new program to estimate the parameters 
of Preston's equation, a general formula for describing the egg shape of birds. \emph{Symmetry} 
15, 231. \doi{10.3390/sym15010231}

Todd, P.H., Smart, I.H.M. (1984) The shape of birds' eggs. \emph{Journal of Theoretical Biology} 
106, 239\eqn{-}243. \doi{10.1016/0022-5193(84)90021-3}
}
\seealso{
\code{\link{EPE}}, \code{\link{lmPE}}, \code{\link{TSE}}
}
\examples{
  zeta <- seq(0, 2*pi, len=2000)
  Par1 <- c(10, 6, 0.325, -0.0415)
  Res1 <- PE(P=Par1, zeta=zeta, simpver=1)
  Par2 <- c(10, 6, -0.325, -0.0415)
  Res2 <- PE(P=Par2, zeta=zeta, simpver=1)

  dev.new()
  plot(Res1$x, Res1$y, asp=1, type="l", col=4, cex.lab=1.5, cex.axis=1.5,
       xlab=expression(italic(x)), ylab=expression(italic(y)))
  lines(Res2$x, Res2$y, col=2)

  dev.new()
  plot(Res1$r, Res2$r, asp=1, cex.lab=1.5, cex.axis=1.5,
       xlab=expression(paste(italic(r), ""[1], sep="")), 
       ylab=expression(paste(italic(r), ""[2], sep="")))
  abline(0, 1, col=4)

  graphics.off()
}