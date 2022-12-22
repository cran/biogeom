\name{GE}
\encoding{UTF-8}
\alias{GE}
\title{
Calculation of the Polar Radius of the Gielis Curve
}
\description{
\code{GE} is used to calculate polar radii of the original Gielis equation  
or one of its simplified versions at given polar angles.
}
\usage{
GE(P, phi, m = 1, simpver = NULL, nval = 1)
}
\arguments{
  \item{P}{the parameters of the original Gielis equation or one of its simplified versions.}
  \item{phi}{the polar angle(s).}
  \item{m}{the given \eqn{m} value that determines the number of angles of 
the Gielis curve within \eqn{[0, 2\pi)}.}
  \item{simpver}{an optional argument to use the simplified version of the original Gielis equation.}
  \item{nval}{the specified value for \eqn{n_{1}} or \eqn{n_{2}} or \eqn{n_{3}} in the simplified versions.}
}
\details{
When \code{simpver = NULL}, the original Gielis equation is selected: 
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{3}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{r} represents the polar radius at the polar angle \eqn{\varphi}; 
\eqn{m} determines the number of angles within \eqn{[0, 2\pi)}; and \eqn{a}, \eqn{k}, \eqn{n_{1}}, 
\eqn{n_{2}}, and \eqn{n_{3}} need to be provided in \code{P}.

\eqn{\quad} When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{a}, \eqn{n_{1}}, and \eqn{n_{2}} need to be provided in \code{P}.

\eqn{\quad} When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}\right)^{-\frac{1}{n_{1}}}, }
where \eqn{a} and \eqn{n_{1}} need to be provided in \code{P}, and \eqn{n_{2}} should be specified in \code{nval}.

\eqn{\quad} When \code{simpver = 3}, the simplified version 3 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}\right)^{-\frac{1}{n_{1}}}, }
where \eqn{a} needs to be provided in \code{P}, and \eqn{n_{1}} 
should be specified in \code{nval}.

\eqn{\quad} When \code{simpver = 4}, the simplified version 4 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}\right)^{-\frac{1}{n_{1}}}, }
where \eqn{a} and \eqn{n_{1}} need to be provided in \code{P}.

\eqn{\quad} When \code{simpver = 5}, the simplified version 5 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{3}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{a}, \eqn{n_{1}}, \eqn{n_{2}}, and \eqn{n_{3}} need to be provided in \code{P}.
 
\eqn{\quad} When \code{simpver = 6}, the simplified version 6 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{a}, \eqn{k}, \eqn{n_{1}}, and \eqn{n_{2}} need to be provided in \code{P}.

\eqn{\quad} When \code{simpver = 7}, the simplified version 7 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{a}, \eqn{k}, and \eqn{n_{1}} need to be provided in \code{P}, and \eqn{n_{2}} 
should be specified in \code{nval}.

\eqn{\quad} When \code{simpver = 8}, the simplified version 8 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{a} and \eqn{k} are parameters that need to be provided in \code{P}, and \eqn{n_{1}} 
should be specified in \code{nval}.

\eqn{\quad} When \code{simpver = 9}, the simplified version 9 is selected:
\deqn{r\left(\varphi\right) = a\left(\left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{1}}\right)^{-\frac{1}{n_{1}}},}
where \eqn{a}, \eqn{k}, and \eqn{n_{1}} need to be provided in \code{P}.
}
\value{
The polar radii predicted by the original Gielis equation or one of its simplified versions.
}
\note{
\code{simpver} here is different from that in the \code{\link{TGE}} function.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Gielis, J. (2003) A generic geometric transformation that unifies a wide range of natural 
and abstract shapes. \emph{American Journal of Botany} 90, 333\eqn{-}338. \doi{10.3732/ajb.90.3.333}

Li, Y., Quinn, B.K., Gielis, J., Li, Y., Shi, P. (2022) Evidence that supertriangles exist 
in nature from the vertical projections of \emph{Koelreuteria paniculata} fruit. \emph{Symmetry} 
14, 23. \doi{10.3390/sym14010023}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Ratkowsky, D.A., Gielis, J. (2020) The generalized Gielis geometric 
equation and its application. \emph{Symmetry} 12, 645. \doi{10.3390/sym12040645}

Shi, P., Xu, Q., Sandhu, H.S., Gielis, J., Ding, Y., Li, H., Dong, X. (2015) Comparison
of dwarf bamboos (\emph{Indocalamus} sp.) leaf parameters to determine relationship
between spatial density of plants and total leaf area per plant. \emph{Ecology and Evolution} 5,
4578\eqn{-}4589. \doi{10.1002/ece3.1728}
}
\seealso{
\code{\link{areaGE}}, \code{\link{curveGE}}, \code{\link{DSGE}}, \code{\link{fitGE}}, 
\code{\link{SurfaceAreaSGE}}, \code{\link{TGE}}, \code{\link{VolumeSGE}}
}
\examples{
GE.par  <- c(2, 1, 4, 6, 3)
varphi.vec <- seq(0, 2*pi, len=2000)
r.theor <- GE(P=GE.par, phi=varphi.vec, m=5)

dev.new()
plot( varphi.vec, r.theor, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic(varphi)), ylab=expression(italic("r")),
      type="l", col=4 ) 

graphics.off()
}