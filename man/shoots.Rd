\name{shoots}
\docType{data}
\alias{shoots}
\title{
Height Growth Data of Bamboo Shoots
}
\description{
The height data of four species of bamboo at Nanjing Forestry University campus in 2016.
}
\usage{
data(shoots)
}
\details{
In the data set, there are four columns of variables: \code{Code}, \code{LatinName}, \code{x}, and \code{y}. 
\code{Code} saves the number codes of different bamboo species; 
\code{LatinName} saves the Latin names of different bamboo species; 
\code{x} saves the investigation times (days from a specific starting time of growth, 
and where every bamboo has a different starting time of growth);  
and \code{y} saves the measured aboveground height values (cm).
 
\code{Code = 1} represents \emph{Phyllostachys iridescens}, and the starting time (namely \code{time = 0}) was defined as 12:00, 3rd April, 2016; 
 
\code{Code = 2} represents \emph{Phyllostachys mannii}, and the starting time (namely \code{time = 0}) was defined as 12:00, 4th April, 2016; 

\code{Code = 3} represents \emph{Pleioblastus maculatus}, and the starting time (namely \code{time = 0}) was defined as 12:00, 29th April, 2016;

\code{Code = 4} represents \emph{Sinobambusa tootsik}, and the starting time (namely \code{time = 0}) was defined as 12:00, 18th April, 2016. 

}
\references{
Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}
}
\examples{
data(shoots)
attach(shoots)
# Choose a species
# 1: Phyllostachys iridescens; 2: Phyllostachys mannii; 
# 3: Pleioblastus maculatus; 4: Sinobambusa tootsik
ind <- 4
x1  <- x[Code == ind]
y1  <- y[Code == ind] 

dev.new()
plot(x1, y1, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="Time (days)", ylab="Height (cm)")

\donttest{
  delta0  <- c(0.5, 1, 2, 5, 10, 20)
  ini.val <- list(600, 25, 0, 40, delta0)
  resu1   <- fitsigmoid(MLRFE, x=x1, y=y1, ini.val=ini.val, simpver=NULL, 
               fig.opt=TRUE, control=list(reltol=1e-20, maxit=20000), 
               subdivisions = 100L, rel.tol=.Machine$double.eps^0.25,                
               abs.tol=.Machine$double.eps^0.25, stop.on.error=TRUE, 
               keep.xy=FALSE, aux=NULL)
}

graphics.off()
}
