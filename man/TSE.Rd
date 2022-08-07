\name{TSE}
\encoding{UTF-8}
\alias{TSE}
\title{
The Todd-Smart Equation (TSE) 
}
\description{
\code{TSE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the Todd and Smart's re-expression of Preston's universal egg shape.
}
\usage{
TSE(P, x, simpver = NULL)
}
\arguments{
  \item{P}{the parameters of the original Todd-Smart equation or one of its simplified versions.}
  \item{x}{the given \eqn{x} values ranging from -1 to 1.}
  \item{simpver}{an optional argument to use the simplified version of the original Todd-Smart equation.}
}
\details{
When \code{simpver = NULL}, the original Preston equation is selected: 
\deqn{y = d_{0}z_{0} + d_{1}z_{1} + d_{2}z_{2} + d_{3}z_{3}}
where 
\deqn{z_{0}=\sqrt{1-x^2}}
\deqn{z_{1}=x\sqrt{1-x^2}}
\deqn{z_{2}=x^{2}\sqrt{1-x^2}}
\deqn{z_{3}=x^{3}\sqrt{1-x^2}}
Here, \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Todd-Smart curve; 
\eqn{d_{0}}, \eqn{d_{1}}, \eqn{d_{2}}, and \eqn{d_{3}} are parameters to be estimated.

\eqn{\quad} When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{y = d_{0}z_{0} + d_{1}z_{1} + d_{2}z_{2}}
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Todd-Smart curve; 
\eqn{d_{0}}, \eqn{d_{1}}, and \eqn{d_{2}} are parameters to be estimated.

\eqn{\quad} When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{y = d_{0}z_{0} + d_{1}z_{1}} 
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Todd-Smart curve; 
\eqn{d_{0}}, and \eqn{d_{1}} are parameters to be estimated.

\eqn{\quad} When \code{simpver = 3}, the simplified version 3 is selected:
\deqn{y = d_{0}z_{0} + d_{2}z_{2}}
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Todd-Smart curve; 
\eqn{d_{0}}, and \eqn{d_{2}} are parameters to be estimated.

}
\value{
  \item{y}{the \eqn{y} value(s) of the Todd-Smart curve corresponding to the given \eqn{x} value(s).}
}
\note{
Here, \eqn{x} and \eqn{y} in the Todd-Smart equation are actually equal to \eqn{y/a} 
and \eqn{x/a}, respectively, in the Preston equation (See \code{\link{PE}} for details). 
Since \eqn{a} represents half the egg length, this means that the egg length is fixed to be 2, 
and the maximum egg width is correspondingly adjusted to keep the same scale.   
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Biggins, J.D., Montgomeries, R.M., Thompson, J.E., Birkhead, T.R. (2022) 
Preston’s universal formula for avian egg shape. \emph{Ornithology} 
In press. \doi{10.1093/ornithology/ukac028}

Biggins, J.D., Thompson, J.E., Birkhead, T.R. (2018) Accurately quantifying 
the shape of birds' eggs. \emph{Ecology and Evolution} 8, 9728\eqn{-}9738. \doi{10.1002/ece3.4412}

Nelder, J.A., Mead, R. (1965). A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Preston, F.W. (1953) The shapes of birds' eggs. \emph{The Auk} 70, 160\eqn{-}182.

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} In press. \doi{10.1111/nyas.14862}

Todd, P.H., Smart, I.H.M. (1984) The shape of birds' eggs. \emph{Journal of Theoretical Biology} 
106, 239\eqn{-}243. \doi{10.1016/0022-5193(84)90021-3}
}
\seealso{
\code{\link{curvePE}}, \code{\link{fitPE}}, \code{\link{lmPE}}, \code{\link{PE}}
}
\examples{
  Par <- c(0.695320398, -0.210538656, -0.070373518, 0.116839895)
  xb1 <- seq(-1, 1, len=20000)
  yb1 <- TSE(P=Par, x=xb1)
  xb2 <- seq(1, -1, len=20000)
  yb2 <- -TSE(P=Par, x=xb2)

  dev.new()
  plot(xb1, yb1, asp=1, type="l", col=2, ylim=c(-1, 1), cex.lab=1.5, cex.axis=1.5, 
    xlab=expression(italic(x)), ylab=expression(italic(y)))
  lines(xb2, yb2, col=4)
}

