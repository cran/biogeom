\name{NRGE}
\encoding{UTF-8}
\alias{NRGE}
\title{
The Narushin-Romanov-Griffin Equation (NRGE)
}
\description{
\code{NRGE} is used to calculate \eqn{y} values at given \eqn{x} values using 
the Narushin-Romanov-Griffin equation (NRGE).
}
\usage{
NRGE(P, x)
}
\arguments{
  \item{P}{the four parameters (i.e., \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D}) of the Narushin-Romanov-Griffin equation.}
  \item{x}{the given \eqn{x} values.}
}
\details{
The Narushin-Romanov-Griffin equation (Narushin et al., 2021) has four parameters in total, 
among which three parameters have clear geometric meanings. 
\deqn{f_{1}(x)=\frac{B}{2}\sqrt{\frac{A^2-4x^2}{A^2+8Cx+4C^2}},}
\deqn{f_{2}(x)=\frac{\sqrt{5.5A^{2}+11AC+4C^{2}} \cdot \left(\sqrt{3}AB-2D\sqrt{A^{2}+2AC+4C^{2
}}\right)}{\sqrt{3}AB\left(\sqrt{5.5A^{2}+11AC+4C^{2}}-2\sqrt{A^{2}+2AC+4C^{2}}\right)},}
\deqn{f_{3}(x)=1-\sqrt{\frac{A\left(A^{2}+8Cx+4C^{2}\right)}{2(A-2C)x^{2}+\left(A^{2}+8AC-4C^{2}\right)x+2AC^{2}+A^{2}C+A^{3}}},}
\deqn{f(x)=\pm f_{1}(x) \cdot \left[1-f_{2}(x) \cdot f_{3}(x)\right].}
Here, \eqn{f(x)} is the Narushin-Romanov-Griffin equation, which is used to predict the \eqn{y} coordinates 
at the given \eqn{x} coordinates; 
\eqn{A} represents the egg length; \eqn{B} represents the egg maximum width; \eqn{C} is a parameter 
to be estimated, and it can be expressed as \eqn{\left(A-B\right)/(2q)}, where \eqn{q} is a parameter to 
be estimated; \eqn{D} represents the egg width associated with \eqn{(3/4)L} from the egg base (to the egg tip) 
on the egg length axis (which can be regarded as the major axis of the egg shape).
}
\value{
The \eqn{y} values predicted by the the Narushin-Romanov-Griffin equation.
}
\note{
Here, parameter \code{C} is a parameter to be estimated, which can be directly calculated numerically based on the egg-shape data.
}
\references{
Narushin, V.G., Romanov, M.N., Griffin, D.K. (2021) Egg and math: introducing a universal formula for egg shape. 
\emph{Annals of the New York Academy of Sciences} 1505, 169\eqn{-}177. \doi{10.1111/nyas.14680}

Shi, P., Gielis, J., Niklas, K.J. (2022) Comparison of a universal (but complex) model for avian egg 
shape with a simpler model. \emph{Annals of the New York Academy of Sciences} 1514, 34\eqn{-}42. \doi{10.1111/nyas.14799}

Tian, F., Wang, Y., Sandhu, H.S., Gielis, J., Shi, P. (2020) Comparison of seed morphology of two ginkgo cultivars. 
\emph{Journal of Forestry Research} 31, 751\eqn{-}758. \doi{10.1007/s11676-018-0770-y}
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\seealso{
\code{\link{curveNRGE}}, \code{\link{fitNRGE}}
}
\examples{

P0 <- c(11.5, 7.8, 1.1, 5.6)
x  <- seq(-11.5/2, 11.5/2, len=2000)
y1 <- NRGE(P=P0, x=x)
y2 <- -NRGE(P=P0, x=x)

dev.new()
plot(x, y1, cex.lab=1.5, cex.axis=1.5, type="l", 
     col=4, ylim=c(-4, 4), asp=1, 
     xlab=expression(italic(x)), ylab=expression(italic(y)) )
lines(x, y2, col=2)

}