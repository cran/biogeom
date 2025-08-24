\name{SAMs}
\encoding{UTF-8}
\docType{data}
\alias{SAMs}
\title{
Boundary Data of Shoot Apical Meristems
}
\description{
The data consist of the boundary data of the shoot apical meristems (SAMs) from six genera.
}
\usage{
data(SAMs)
}
\details{
In the data set, there are three columns of variables: \code{Genus}, \code{x}, and \code{y}. 
\code{Genus} saves the genus information for each SAM; 
\code{x} saves the \eqn{x} coordinates of the six SAMs in the Cartesian coordinate system (unitless); 
and \code{y} saves the \eqn{y} coordinates of the six SAMs in the Cartesian coordinate system (unitless). 
In \code{Genus}, there are six genera including \emph{Glycine}, \emph{Coleus}, \emph{Juniperus}, \emph{Pinus}, 
\emph{Equisetum}, and \emph{Isoetes}, spanning the three lineages: angiosperms, gymnosperms, and pteridophytes.
}
\references{
Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}
}
\seealso{
\code{\link{fitAM}}, \code{\link{PlanCoor}}, \code{\link{SurfaceAreaAM}}, \code{\link{VolumeAM}}
}
\examples{
#### See Shi et al. (2025) for details #########################################
# Shi, P., Liu, X., Gielis, J., Beirinckx, B., Niklas, K.J. (2025) 
#     Comparison of six non-linear equations in describing the 2-D 
#     profiles of apical meristems. American Journal of Botany (under review).
################################################################################

data(SAMs)

uni.sam <- sort( unique(SAMs$Genus) )
ind     <- 2
Data    <- SAMs[SAMs$Genus==uni.sam[ind], ]
x0      <- Data$x
y0      <- Data$y

Res1    <- adjdata(x0, y0, ub.np=200, times=1.2, len.pro=1/20)
X       <- Res1$x
Y       <- Res1$y

dev.new()
plot( X, Y, pch=1, cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(paste(italic(x), " (unitless)", sep="")), 
      ylab=expression(paste(italic(y), " (unitless)", sep="")) )

graphics.off()
}
