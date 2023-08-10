\name{DNRGE}
\encoding{UTF-8}
\alias{DNRGE}
\title{
Calculation of the First-Order Derivative of the Narushin-Romanov-Griffin Equation  
}
\description{
\code{DNRGE} is used to calculate the first-order derivative of the Narushin-Romanov-Griffin equation at a given \emph{x}-value.
}
\usage{
DNRGE(P, x)
}
\arguments{
  \item{P}{the parameters of the Narushin-Romanov-Griffin equation.}
  \item{x}{the \emph{x}-value used in the Narushin-Romanov-Griffin equation.}
}
\details{
Let us define:

\deqn{f_{1}(x) = \frac{B}{2}\sqrt{\frac{A^2-4x^2}{A^2+8Cx+4C^2}},}
\deqn{f_{2}(x) = \sqrt{\frac{A\left(A^{2}+8Cx+4C^{2}\right)}{2(A-2C)x^{2}+\left(A^{2}+8AC-4C^{2}\right)x+2AC^{2}+A^{2}C+A^{3}}},}
\deqn{f_{3}(x) = A^2 - 4x,}
\deqn{f_{4}(x) = A^2+8Cx+4C^2,}
\deqn{E = \frac{\sqrt{5.5A^{2}+11AC+4C^{2}} \cdot \left(\sqrt{3}AB-2D\sqrt{A^{2}+2AC+4C^{2
}}\right)}{\sqrt{3}AB\left(\sqrt{5.5A^{2}+11AC+4C^{2}}-2\sqrt{A^{2}+2AC+4C^{2}}\right)},}
\deqn{F = 2\left(A-2C\right),}
\deqn{G = A^{2}+8AC-4C^{2},}
\deqn{H = 2AC^{2}+A^{2}C+A^{3},}

and then the first-order derivative of the Narushin-Romanov-Griffin equation at a given \emph{x}-value is:
\deqn{J(x) = -\frac{4\,f_{1}(x)\left[C\,f_{3}(x)+x\,f_{4}(x)\right]}{f_{3}(x) \cdot f_{4}(x)}\left\{1-E \cdot \left[1-f_{2}(x)\right]\right\}-\frac{AE}{2}\frac{f_{1}(x)}{f_{2}(x)}\frac{f_{4}(x) \cdot \left(2Fx+G\right)}{\left(Fx^2+Gx+H\right)^2},}
where \code{P} has four parameters: \eqn{A}, \eqn{B}, \eqn{C}, and \eqn{D}.
}
\note{
The argument \code{P} in the \code{\link{DNRGE}} function has the same parameters, as those in the 
\code{\link{NRGE}} function.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Narushin, V.G., Romanov, M.N., Griffin, D.K. (2021) Egg and math: introducing a universal formula for egg shape. 
\emph{Annals of the New York Academy of Sciences} 1505, 169\eqn{-}177. \doi{10.1111/nyas.14680}

Narushin, V.G., Romanov, M.N., Mishra, B., Griffin, D.K. (2022) Mathematical progression of 
avian egg shape with associated area and volume determinations. 
\emph{Annals of the New York Academy of Sciences} 1513, 65\eqn{-}78. \doi{10.1111/nyas.14771}

Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2023) 
A simple way to calculate the volume and surface area of avian eggs. 
\emph{Annals of the New York Academy of Sciences} 1524, 118\eqn{-}131. \doi{10.1111/nyas.15000}
}
\seealso{
\code{\link{fitNRGE}}, \code{\link{NRGE}}, \code{\link{SurfaceAreaNRGE}}
}
\examples{
  Par6 <- c(4.51, 3.18, 0.1227, 2.2284)
  xx3  <- seq(-4.51/2, 4.51/2, len=2000)
  J1   <- DNRGE(P=Par6, x=xx3)
  J2   <- -DNRGE(P=Par6, x=xx3)
  ind  <- which(is.na(J1) | is.na(J2))
  xx3  <- xx3[-ind]
  J1   <- J1[-ind]
  J2   <- J2[-ind]

  dev.new()
  plot(xx3, J1, type="l", col=4, cex.lab=1.5, cex.axis=1.5,
       xlim=c(-4.51/2, 4.51/2), ylim=c(-20, 20), xlab=expression(italic(x)), 
       ylab=expression(paste(italic(J), "(", italic(x), ")", sep="")))
  lines(xx3, J2, col=2)  

  graphics.off()
}