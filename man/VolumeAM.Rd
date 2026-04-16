\name{VolumeAM}
\encoding{UTF-8}
\alias{VolumeAM}
\title{
Calculation of the Volume of an Apical Meristem. 
}
\description{
\code{VolumeAM} is used to calculate the volume of an apical meristem that 
follows the the hybrid catenary-parabolic equation or the superparabolic equation.
}
\usage{
VolumeAM(model = "Hybrid", P, upper = Inf, subdivisions = 100L,
         rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
         stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{model}{an optional argument to select a model to calculate the volume.}
  \item{P}{the parameters of the hybrid catenary-parabolic equation or the superparabolic equation.}
  \item{upper}{the upper limit of integration in \eqn{x}.}
  \item{subdivisions}{please see the arguments for the \code{\link[stats]{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link[stats]{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link[stats]{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link[stats]{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link[stats]{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link[stats]{integrate}} function in package \pkg{stats}.}
}
\details{
The formula of the volume (\eqn{V}) of an apical meristem 
based on the hybrid catenary-parabolic equation or the superparabolic equation have two cases. 

\eqn{\quad} Case (i): if the AM's profile is a downward-opening curve, \eqn{V} takes the form 
\deqn{V(x) = 2\, \pi \int_{0}^{a} x\, \left|y\left(x\right) \right| dx,} 
where \eqn{a} is the upper limit of integration in \eqn{x}. The lower limit of integration in \eqn{x} is 0. 

\eqn{\quad} Case (ii): if the AM's profile is an upward-opening curve, \eqn{V} takes the form 
\deqn{V(x) = \pi a^2 \left|y\left(a\right)\right| - 2\, \pi \int_{0}^{a} x\, \left|y\left(x\right) \right| dx,} 
where \eqn{a} is the upper limit of integration in \eqn{x}. The lower limit of integration in \eqn{x} is 0. 

\eqn{\quad} If \code{model = "Hybrid"}, \eqn{y} denotes the hybrid catenary-parabolic equation, 
which equals \eqn{y\left(x\right) = \alpha\,\mbox{cosh}(\beta x) + \gamma x^{2}-\alpha}. 
Here, \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma} are model parameters provided by the argument \code{P}.

\eqn{\quad} If \code{model = "Superparabola"}, \eqn{y} denotes the superparabolic equation, 
which equals \eqn{y\left(x\right) = \beta_{1}\, {\left|x\right|}^{\beta_{2}}}. 
Here, \eqn{\beta_{1}} and \eqn{\beta_{2}} are model parameters provided by the argument \code{P}.  
}
\note{
The argument \code{P} in the \code{\link{VolumeAM}} function has three parameters for \code{model = "Hybrid"}, 
and two parameters for \code{model = "Superparabola"}. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Liu, X., Gielis, J., Beirinckx, B., Niklas, K.J. (2026) 
Comparison of six nonlinear equations describing the 2-D profiles of apical meristems. 
\emph{American Journal of Botany} 113, e70177. \doi{10.1002/ajb2.70177}
}
\seealso{
\code{\link{fitAM}}, \code{\link{PlanCoor}}, \code{\link{SurfaceAreaAM}}, \code{\link{SAMs}}
}
\examples{
  par3 <- c(-20, 0.9695, 11.40)
  VolumeAM(model = "Hybrid", P = par3, upper=1)
  
  par5 <- c(1.3807, 1.3173)
  VolumeAM(model = "Superparabola", P = par5, upper=1)
}

