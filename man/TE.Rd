\name{TE}
\encoding{UTF-8}
\alias{TE}
\title{
The Troscianko Equation (TE) 
}
\description{
\code{TE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the re-expression of Troscianko's egg-shape equation, which was proposed by Biggins et al. (2018, 2022).
}
\usage{
TE(P, x)
}
\arguments{
  \item{P}{the parameters of the Troscianko equation, including \eqn{\alpha_{0}}, \eqn{\alpha_{1}}, and \eqn{\alpha_{2}}.}
  \item{x}{the given \eqn{x} values ranging from \eqn{-}1 to 1.}
}
\details{
The Troscianko equation is recommended as (Biggins et al., 2022): 
\deqn{y = \exp\left(\alpha_{0}+\alpha_{1}\,x+\alpha_{2}\,x^2\right)\sqrt{1-x^2}, }
where \eqn{x} and \eqn{y} represent the abscissa and ordinate of an arbitrary point on the Troscianko curve; 
\eqn{\alpha_{0}}, \eqn{\alpha_{1}}, and \eqn{\alpha_{2}} are parameters to be estimated.
}
\value{
The \eqn{y} values predicted by the Troscianko equation.
}
\note{
Here, \eqn{x} and \eqn{y} in the Troscianko equation are actually equal to \eqn{y/a} 
and \eqn{x/a}, respectively, in the explicit Troscianko equation, where \eqn{a} 
represents half the egg length (See \code{\link{ETE}} for details). 
This means that the egg length is scaled to be 2, 
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

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Wang, L., Quinn, B.K., Gielis, J. (2022) A new program to estimate the parameters 
of Preston's equation, a general formula for describing the egg shape of birds. \emph{Symmetry} 
Under review.

Troscianko, J. (2014). A simple tool for calculating egg shape, volume and surface area from digital images. 
\emph{Ibis}, 156, 874\eqn{-}878. \doi{10.1111/ibi.12177} 
}
\seealso{
\code{\link{fitETE}}, \code{\link{lmTE}}
}
\examples{
  Par <- c(-0.377, -0.29, -0.16)
  xb1 <- seq(-1, 1, len=20000)
  yb1 <- TE(P=Par, x=xb1)
  xb2 <- seq(1, -1, len=20000)
  yb2 <- -TE(P=Par, x=xb2)

  dev.new()
  plot(xb1, yb1, asp=1, type="l", col=2, ylim=c(-1, 1), cex.lab=1.5, cex.axis=1.5, 
    xlab=expression(italic(x)), ylab=expression(italic(y)))
  lines(xb2, yb2, col=4)

  graphics.off()
}

