\name{LeafSizeDist}
\encoding{UTF-8}
\docType{data}
\alias{LeafSizeDist}
\title{
Leaf size distribution of \emph{Shibataea chinensis}
}
\description{
The data consist of the leaf size measures of 12 individuad culms of \emph{Shibataea chinensis}.
}
\usage{
data(LeafSizeDist)
}
\details{
In the data set, there are five columns of variables: \code{Code}, \code{L}, 
\code{W}, \code{A}, and \code{M}. \code{Code} saves the bamboo number (the number before the hyphen) 
and the leaf number (the number after the hyphen) on each bamboo. \code{L} saves the length of each 
leaf lamina (cm); \code{W} saves the width of each leaf lamina (cm); 
\code{A} saves the area of each leaf lamina (\eqn{\mbox{cm}^{2}}); and \code{M} saves 
the dry mass of each leaf lamina (g).
}
\references{
Lian, M., Shi, P., Zhang, L., Yao, W., Gielis, J., Niklas, K.J. (2023) A generalized performance equation 
and its application in measuring the Gini index of leaf size inequality. 
\emph{Trees \eqn{-} Structure and Function} 37, 1555\eqn{-}1565. \doi{10.1007/s00468-023-02448-8}
}
\examples{
data(LeafSizeDist)

CulmNumber <- c()
for(i in 1:length(LeafSizeDist$Code)){
  temp <- as.numeric( strsplit(LeafSizeDist$Code[i],"-", fixed=TRUE)[[1]][1] )
  CulmNumber <- c(CulmNumber, temp)
}

uni.CN <- sort( unique(CulmNumber) )  
ind    <- CulmNumber==uni.CN[1]

A0 <- LeafSizeDist$A[ind]
A1 <- sort( A0 )
x1 <- 1:length(A1)/length(A1)
y1 <- cumsum(A1)/sum(A1)
x1 <- c(0, x1)
y1 <- c(0, y1)

M0 <- LeafSizeDist$M[ind]
M1 <- sort( M0 )
x2 <- 1:length(M1)/length(M1)
y2 <- cumsum(M1)/sum(M1)
x2 <- c(0, x2)
y2 <- c(0, y2)

dev.new()
plot(x1, y1, cex.lab=1.5, cex.axis=1.5, type="l",
    xlab="Cumulative proportion of the number of leaves", 
    ylab="Cumulative proportion of leaf area") 
lines(c(0,1), c(0,1), col=4) 

dev.new()
plot(x2, y2, cex.lab=1.5, cex.axis=1.5, type="l",
    xlab="Cumulative proportion of the number of leaves", 
    ylab="Cumulative proportion of leaf dry mass")
lines(c(0,1), c(0,1), col=4)   

graphics.off()
}
