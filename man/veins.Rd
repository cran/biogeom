\name{veins}
\docType{data}
\encoding{UTF-8}
\alias{veins}
\title{
Leaf Vein Data of \emph{Michelia compressa}
}
\description{
The data consist of the leaf vein data of a leaf of \emph{M. compressa} 
sampled at Nanjing Forestry University campus in late July 2019.
}
\usage{
data(veins)
}
\details{
In the data set, there are two columns of variables: \code{x} and \code{y}. 
\code{x} saves the \eqn{x} coordinates of the leaf veins in the Cartesian coordinate system (cm); 
\code{y} saves the \eqn{y} coordinates of the leaf veins in the Cartesian coordinate system (cm).
}
\note{
The data cannot be adjusted by the \code{\link{adjdata}} function.
}
\seealso{
\code{\link{fracdim}}
}
\references{
Shi, P., Yu, K., Niinemets, Ü., Gielis, J. (2021) Can leaf shape be represented by the ratio of 
leaf width to length? Evidence from nine species of \emph{Magnolia} and \emph{Michelia} (Magnoliaceae). 
\emph{Forests} 12, 41. \doi{10.3390/f12010041}
}
\examples{
data(veins)

dev.new()
plot(veins$x, veins$y, cex=0.01, asp=1, cex.lab=1.5, cex.axis=1.5, 
     xlab=expression(italic("x")), ylab=expression(italic("y")))

}
