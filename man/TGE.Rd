\name{TGE}
\encoding{UTF-8}
\alias{TGE}
\title{
Calculation of the Polar Radius of the Twin Gielis Curve
}
\description{
\code{TGE} is used to calculate the polar radii of the twin Gielis equation or 
one of its simplified versions at given polar angles.
}
\usage{
TGE(P, phi, m = 1, simpver = NULL, nval = 1)
}
\arguments{
  \item{P}{the parameters of the twin Gielis equation or one of its simplified versions.}
  \item{phi}{the polar angle(s).}
  \item{m}{the given \eqn{m} value that determines the number of angles of the twin Gielis curve within \eqn{[0, 2\pi)}.}
  \item{simpver}{an optional argument to use the simplified version of the twin Gielis equation.}
  \item{nval}{the specified value for \eqn{n_{2}} or \eqn{n_{3}} in the simplified versions.}
}
\details{
The general form of the twin Gielis equation can be represented as follows:

\deqn{r\left(\varphi\right) = \mathrm{exp}\left\{\frac{1}{\alpha+\beta\,\mathrm{ln}\left[r_{e}\left(\varphi\right)\right]}+\gamma\right\},}
where \eqn{r} represents the polar radius of the twin Gielis curve at the polar angle \eqn{\varphi}, and 
\eqn{r_{e}} represents the elementary polar radius at the polar angle \eqn{\varphi}. There is a hyperbolic 
link function to link their log-transformations, i.e., 
\deqn{\mathrm{ln}\left[r\left(\varphi\right)\right] = \frac{1}{\alpha+\beta\,\mathrm{ln}\left[r_{e}\left(\varphi\right)\right]}+\gamma.}  
The first three elements of \code{P} are \eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma}, and the remaining element(s) of 
\code{P} are the parameters of the elementary polar function, i.e., \eqn{r_{e}\left(\varphi\right)}. 
See Shi et al. (2020) for details. 

\eqn{\quad} When \code{simpver = NULL}, the original twin Gielis equation is selected: 
\deqn{r_{e}\left(\varphi\right) = \left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{3}},}
where \eqn{r_{e}} represents the elementary polar radius at the polar angle \eqn{\varphi}; 
\eqn{m} determines the number of angles of the twin Gielis curve within \eqn{[0, 2\pi)}; 
and \eqn{k}, \eqn{n_{2}}, and \eqn{n_{3}} are the fourth to the sixth elements in \code{P}. In total, there are 
six elements in \code{P}.

\eqn{\quad} When \code{simpver = 1}, the simplified version 1 is selected:
\deqn{r_{e}\left(\varphi\right) = \left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}},}
where \eqn{n_{2}} is the fourth element in \code{P}. There are 
four elements in total in \code{P}.

\eqn{\quad} When \code{simpver = 2}, the simplified version 2 is selected:
\deqn{r_{e}\left(\varphi\right) = \left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}, }
where \eqn{n_{2}} should be specified in \code{nval}, and \code{P} only includes three elements, i.e., 
\eqn{\alpha}, \eqn{\beta}, and \eqn{\gamma}.

\eqn{\quad} When \code{simpver = 3}, the simplified version 3 is selected:
\deqn{r_{e}\left(\varphi\right) = \left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{3}},}
where \eqn{n_{2}} and \eqn{n_{3}} are the fourth and fifth elements in \code{P}. There are 
five elements in total in \code{P}.
 
\eqn{\quad} When \code{simpver = 4}, the simplified version 4 is selected:
\deqn{r_{e}\left(\varphi\right) = \left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}},}
where \eqn{k} and \eqn{n_{2}} are the fourth and fifth elelments in \code{P}. There are 
five elements in total in \code{P}.
 
\eqn{\quad} When \code{simpver = 5}, the simplified version 5 is selected:
\deqn{r_{e}\left(\varphi\right) = \left|\mathrm{cos}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}}+
\left|\frac{1}{k}\,\mathrm{sin}\left(\frac{m}{4}\varphi\right)\right|^{n_{2}},}
where \eqn{k} is the fourth elelment in \code{P}. There are 
four elements in total in \code{P}. \eqn{n_{2}} should be specified in \code{nval}.

}
\value{
The polar radii predicted by the twin Gielis equation or one of its simplified versions.
}
\note{
\code{simpver} here is different from that in the \code{\link{GE}} function.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Li, Y., Quinn, B.K., Gielis, J., Li, Y., Shi, P. (2022) Evidence that supertriangles exist 
in nature from the vertical projections of \emph{Koelreuteria paniculata} fruit. \emph{Symmetry} 
14, 23. \doi{10.3390/sym14010023}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Ratkowsky, D.A., Gielis, J. (2020) The generalized Gielis geometric 
equation and its application. \emph{Symmetry} 12, 645. \doi{10.3390/sym12040645}
}
\seealso{
\code{\link{areaGE}}, \code{\link{curveGE}}, \code{\link{fitGE}}, \code{\link{GE}}
}
\examples{
TGE.par    <- c(2.88, 0.65, 1.16, 139)
varphi.vec <- seq(0, 2*pi, len=2000)
r2.theor   <- TGE(P=TGE.par, phi=varphi.vec, simpver=1, m=5)

dev.new()
plot( varphi.vec, r2.theor, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic(varphi)), ylab=expression(italic("r")),
      type="l", col=4 ) 

starfish4 <- curveGE(TGE, P=c(0, 0, 0, TGE.par), simpver=1, m=5, fig.opt=TRUE)
}