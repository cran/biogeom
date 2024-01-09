curveovate <- function(expr, P, x, fig.opt = FALSE, 
                deform.fun = NULL, Par = NULL,
                xlim = NULL, ylim = NULL, unit = NULL, main = NULL){
  if((is.null(deform.fun) & !is.null(Par
      )) | (!is.null(deform.fun) & is.null(Par)))
  stop("'Par' should be provided when 'deform.fun' is not null.")
  p        <- length(P)   
  x0       <- P[1]
  y0       <- P[2]
  theta    <- P[3]
  PAR      <- P[4:p] 
  UpperFun <- function(P, x){
    z1    <- sort(x, decreasing=TRUE, index.return=TRUE)$x
    temp1 <- expr(P, z1, simpver=1) 
    Temp  <- cbind(z1, temp1)
    Temp  <- na.omit(Temp)
    x     <- Temp[,1]
    y     <- Temp[,2]
    list(x=x, y=y)
  }
  LowerFun <- function(P, x){
    z1    <- sort(x, decreasing=FALSE, index.return=TRUE)$x
    temp1 <- expr(P, z1, simpver=1) 
    temp2 <- -temp1
    Temp  <- cbind(z1, temp2)
    Temp  <- na.omit(Temp)
    x     <- Temp[,1]
    y     <- Temp[,2]
    list(x=x, y=y)
  }
  res1     <- UpperFun(PAR, x)
  res2     <- LowerFun(PAR, x)
  x1       <- c(res1$x, res2$x)
  y1       <- c(res1$y, res2$y)
  if(!is.null(deform.fun)){  
    Resu <- deform.fun(Par=Par, z=cbind(x1,y1))
    x1   <- Resu$x
    y1   <- Resu$y
  }
  x.rot    <- x1*cos(theta) - y1*sin(theta)
  y.rot    <- y1*cos(theta) + x1*sin(theta)   
  x.coordi <- x.rot + x0
  y.coordi <- y.rot + y0    
  if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True"){
        if(is.null(xlim)) xlim <- NULL
        if(is.null(ylim)) ylim <- NULL
        if(!is.null(xlim)) xlim <- xlim
        if(!is.null(ylim)) ylim <- ylim
        if(!is.null(unit)){
  xlabel <- bquote(paste(italic("x")," (", .(unit),")", sep="")) 
  ylabel <- bquote(paste(italic("y")," (", .(unit),")", sep=""))
        }
        if(is.null(unit)){
          xlabel <- bquote( bolditalic("x") ) 
          ylabel <- bquote( bolditalic("y") )
        }
        dev.new()
        plot( x.coordi, y.coordi, asp=1, xlab=xlabel, ylab=ylabel, 
              pch=1, cex=2, cex.lab=1.5,
              cex.axis=1.5, type="l", lwd=2 )
        abline(h=y0, lty=2, col=4)
        abline(v=x0, lty=2, col=4)   
        abline(y0-tan(theta)*x0, tan(theta), col=2, lty=2)
        title(main=main, cex.main=1.5, col.main=4, font.main=1)
  }
  list(x = x.coordi, y = y.coordi)
}
