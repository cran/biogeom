\name{fracdim}
\encoding{UTF-8}
\alias{fracdim}
\title{
Calculation of Fractal Dimension of Lef Veins Based on the Box-Counting Method
}
\description{
\code{fracdim} is used to calculate the fractal dimension of leaf veins 
based on the box-counting method.
}
\usage{
fracdim(x, y, frac.fig = TRUE, denomi.range = seq(8, 30, by=1), 
        ratiox = 0.02, ratioy = 0.08, main = NULL)
}
\arguments{
  \item{x}{the \eqn{x} coordinates of leaf-vein pixels.}
  \item{y}{the \eqn{y} coordinates of leaf-vein pixels.}
  \item{frac.fig}{the option of drawing the results of the linear fitting.}
  \item{denomi.range}{the number of equidistant segments of the maximum range 
between the range of the \eqn{x} coordinates and that of the \eqn{y} coordinates.}
  \item{ratiox}{the the \eqn{x} coordinate of the location parameter for positioning the legend.}
  \item{ratioy}{the the \eqn{y} coordinate of the location parameter for positioning the legend.}
  \item{main}{the main title of the figure.}
}
\details{
The box-counting approach uses a group of boxes (squares for simplicity) with different
sizes (\eqn{\delta}) to divide the leaf vein image into different parts. Let \eqn{N} represent the number
of boxes that include at least one pixel of leaf vein. 
The maximum of the range of the \eqn{x} coordinates and the range of the \eqn{y} coordinates 
for leaf-vein pixels is defined as \eqn{z}. Let \eqn{\delta} represent the vector of 
\eqn{z}/\code{denomi.range}. Then, we used the following equation to calculate the fractal 
dimension of leaf veins:
\deqn{\mathrm{ln } N = a + b \mathrm{ ln} {\delta}^{-1},}
where \eqn{b} is the theoretical value of the fractal dimension. We can use its estimate as the 
numerical value of the fractal dimension for a leaf venation network.      
}
\value{
  \item{a}{the estimate of the intercept.}
  \item{sd.a}{the standard deviation of the estimated intercept.}
  \item{lci.a}{the lower bound of the 95\% confidence interval of the estimated intercept.}
  \item{uci.a}{the upper bound of the 95\% confidence interval of the estimated intercept.}
  \item{b}{the estimate of the slope.}
  \item{sd.b}{the standard deviation of the estimated slope.}
  \item{lci.a}{the lower bound of the 95\% confidence interval of the estimated slope.}
  \item{uci.a}{the upper bound of the 95\% confidence interval of the estimated slope.}
  \item{r.sq}{the coefficient of determination.}
  \item{delta}{the vector of box sizes.}
  \item{N}{the number of boxes that include at least one pixel of leaf vein.}
}
\note{
Here, \code{x} and \code{y} cannot be adjusted by the \code{\link{adjdata}} function 
because the leaf veins are not the leaf's boundary data.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Yu, K., Niinemets, Ü., Gielis, J. (2021) Can leaf shape be represented by the ratio of 
leaf width to length? Evidence from nine species of \emph{Magnolia} and \emph{Michelia} (Magnoliaceae). 
\emph{Forests} 12, 41. \doi{10.3390/f12010041}

Vico, P.G., Kyriacos, S., Heymans, O., Louryan, S., Cartilier, L. (1998) 
Dynamic study of the extraembryonic vascular network of the 
chick embryo by fractal analysis. \emph{Journal of Theoretical Biology} 195, 525\eqn{-}532. 
\doi{10.1006/jtbi.1998.0810}
}
\seealso{
\code{\link{veins}}
}
\examples{
data(veins)

dev.new()
plot(veins$x, veins$y, cex=0.01, asp=1, cex.lab=1.5, cex.axis=1.5, 
     xlab=expression(italic("x")), ylab=expression(italic("y")))

fracdim(veins$x, veins$y)
}