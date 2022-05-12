curveNRGE <- function(P, x, fig.opt = FALSE, deform.fun = NULL, 
    Par = NULL, xlim = NULL, ylim = NULL, unit = NULL, main = NULL){

  if( length(P)!=7) stop("There are seven elements in P!" )  
  if((is.null(deform.fun) & !is.null(Par)) | (!is.null(deform.fun) & is.null(Par)))
      stop("'Par' should be provided when 'expr' is not null.")

  UpperFun <- function(P, x){
    A    <- P[1]
    B    <- P[2]
    C    <- P[3]
    D    <- P[4]
    fun0 <- function(z){
        z[z < -A/2] <- -A/2
        z[z > A/2]  <- A/2
        return(z)
    }
    x     <- fun0(x)
    z1    <- sort(x, decreasing=TRUE, index.return=TRUE)$x
    temp1 <- B/2*sqrt((A^2-4*z1^2)/(A^2+8*C*z1+4*C^2))*
             (1-(sqrt(5.5*A^2+11*B*C+4*C^2)*(
             sqrt(3)*B*A-2*D*sqrt(A^2+2*C*A+4*C^2)))/
             (sqrt(3)*B*A*(sqrt(5.5*A^2+11*A*C+4*C^2
             )-2*sqrt(A^2+2*C*A+4*C^2)))*
             (1-sqrt((A*(A^2+8*C*z1+4*C^2))/(2*(A-2*C)*z1^2+
             (A^2+8*A*C-4*C^2)*z1+2*A*C^2+A^2*C+A^3)))) 

    Temp  <- cbind(z1, temp1)
    Temp  <- na.omit(Temp)
    list(x=Temp[,1], y=Temp[,2])
  }


  LowerFun <- function(P, x){
    A     <- P[1]
    B     <- P[2]
    C     <- P[3]
    D     <- P[4]
    fun0 <- function(z){
        z[z < -A/2] <- -A/2
        z[z > A/2]  <- A/2
        return(z)
    }
    x     <- fun0(x)
    z1    <- sort(x, decreasing=FALSE, index.return=TRUE)$x
    temp1 <- B/2*sqrt((A^2-4*z1^2)/(A^2+8*C*z1+4*C^2))*
             (1-(sqrt(5.5*A^2+11*B*C+4*C^2)*(sqrt(3
             )*B*A-2*D*sqrt(A^2+2*C*A+4*C^2)))/
             (sqrt(3)*B*A*(sqrt(5.5*A^2+11*A*C+4*C^2
             )-2*sqrt(A^2+2*C*A+4*C^2)))*
             (1-sqrt((A*(A^2+8*C*z1+4*C^2))/(2*(A-2*C)*z1^2+
             (A^2+8*A*C-4*C^2)*z1+2*A*C^2+A^2*C+A^3))))    
    temp2 <- -temp1
    Temp  <- cbind(z1, temp2)
    Temp  <- na.omit(Temp)
    x     <- Temp[,1]
    y     <- Temp[,2]
    list(x=x, y=y)
  }


  x0       <- P[1]
  y0       <- P[2]
  theta    <- P[3]
  A        <- P[4]
  B        <- P[5]
  C        <- P[6]
  D        <- P[7]
  PAR      <- c(A, B, C, D)  
  res1     <- UpperFun(P=PAR, x=x)
  res2     <- LowerFun(P=PAR, x=x)
  x        <- c(res1$x, res2$x)
  y        <- c(res1$y, res2$y)

  if(!is.null(deform.fun)){  
    Resu <- deform.fun(Par=Par, z=cbind(x,y))
    x    <- Resu$x
    y    <- Resu$y
  }

  x.rot    <- x*cos(theta) - y*sin(theta)
  y.rot    <- y*cos(theta) + x*sin(theta) 
  
  x.coordi <- x.rot + x0
  y.coordi <- y.rot + y0


  res3 <- UpperFun(P=PAR, x=A)
  xu <- res3$x
  yu <- res3$y

  if(!is.null(deform.fun)){  
    Right <- deform.fun(Par=Par, z=cbind(xu,yu))
    xu    <- Right$x
    yu    <- Right$y
  }    
    
  xv <- xu*cos(theta) - yu*sin(theta)
  yv <- yu*cos(theta) + xu*sin(theta) 
  xv <- xv + x0
  yv <- yv + y0 
    

  if(fig.opt == "T" | fig.opt == "TRUE"){
        if(is.null(xlim)) xlim <- NULL
        if(is.null(ylim)) ylim <- NULL
        if(!is.null(xlim)) xlim <- xlim
        if(!is.null(ylim)) ylim <- ylim
        if(!is.null(unit)){
    xlabel <- bquote(paste(italic("x")," (", .(unit), ")",sep="")) 
    ylabel <- bquote(paste(italic("y")," (", .(unit), ")",sep=""))
        }
        if(is.null(unit)){
          xlabel <- bquote( italic("x") ) 
          ylabel <- bquote( italic("y") )
        }
        dev.new()
        plot( x.coordi, y.coordi, asp=1, xlab=xlabel, ylab=ylabel, 
              pch=1, cex=1, cex.lab=1.5, cex.axis=1.5, lwd=2, 
              type="l", xlim=xlim, ylim=ylim)
        abline(h=y0, lty=2, col=4)
        abline(v=x0, lty=2, col=4)   
        abline(y0-tan(theta)*x0, tan(theta), col=2, lty=2)
        points(x0, y0, cex=2, pch=16, col=2)
        points(xv, yv, cex=2, pch=1, col=4)
        title(main=main, cex.main=1.5, col.main=4, font.main=1)
  }
  list(x = x.coordi, y = y.coordi)
}