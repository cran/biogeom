\name{sigmoid}
\encoding{UTF-8}
\alias{sigmoid}
\title{
Sigmoid Growth Equation
}
\description{
\code{sigmoid} is used to calculate the \eqn{y} values (e.g., biomass, height, body length, and so on) at given investigation times.
}
\usage{
sigmoid(expr, P, x, simpver = 1, subdivisions = 100L,
        rel.tol = .Machine$double.eps^0.25, 
        abs.tol = rel.tol, stop.on.error = TRUE, 
        keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{expr}{a performance equation or one of its simplified versions.}
  \item{P}{the parameters of the performance equation or one of its simplified versions.}
  \item{x}{the given investigation times.}
  \item{simpver}{an optional argument to use the simplfied version of the performance equation.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
}
\details{
The performance equations denote \code{MbetaE}, \code{MBriereE}, 
\code{MLRFE}, and their simplified versions. 
The arguments of \code{P} and \code{simpver} should correspond 
to \code{expr} (i.e., \code{\link{MbetaE}}, \code{\link{MBriereE}}, 
\code{\link{MLRFE}}, and \code{\link{MPerformanceE}}). 
The sigmoid curve is the integral of the performance equation or one of its simplified versions.
}
\value{
The \eqn{y} values (i.e., biomass, height, body length, and so on) at given investigation times.  
The growth euqation is actually an integral of the performance equation or one of its simplified versions.
}
\note{
Here, the user can define other performance equations, but new equations or 
their simplified versions should include the lower and upper thresholds in 
the \eqn{x}-axis corresponding to \eqn{y = 0}, whose indices of the parameters in \code{P} should 
be the same as those in \code{\link{MbetaE}}, \code{\link{MBriereE}}, 
\code{\link{MLRFE}}, and \code{\link{MPerformanceE}}).
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Jin, J., Quinn, B.K., Shi, P. (2022) The modified Brière equation and its 
applications. \emph{Plants} 11, 1769. \doi{10.3390/plants11131769}

Lian, M., Shi, P., Zhang, L., Yao, W., Gielis, J., Niklas, K.J. (2023) A generalized performance equation 
and its application in measuring the Gini index of leaf size inequality. 
\emph{Trees \eqn{-} Structure and Function} 37, 1555\eqn{-}1565. \doi{10.1007/s00468-023-02448-8}

Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}
}
\seealso{
\code{\link{fitsigmoid}}, \code{\link{MbetaE}}, \code{\link{MBriereE}}, 
\code{\link{MLRFE}}, and \code{\link{MPerformanceE}}
}
\examples{

Pa1 <- c(3, 3, 10, 1)
xv1 <- seq(-5, 15, len=2000)
yv1 <- sigmoid(MBriereE, P=Pa1, x=xv1, simpver=1)
Pa2 <- c(3, 3, 2, 12, 1)
yv2 <- sigmoid(MBriereE, P=Pa2, x=xv1, simpver=NULL)
Pa3 <- c(200, 1.32, 0.85, 12, 3.7, 17)
yv3 <- sigmoid(MPerformanceE, P=Pa3, x=xv1, simpver=1)

dev.new()
plot( xv1, yv2, cex.lab=1.5, cex.axis=1.5, type="l", col=4,
      xlab=expression(italic(x)), ylab=expression(italic(y)) )
lines( xv1, yv1, col=2 )
lines( xv1, yv3, col=3 )

graphics.off()
}