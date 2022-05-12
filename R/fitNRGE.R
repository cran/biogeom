

fitNRGE <- function(x, y, angle = NULL, x0=NULL, y0=NULL,
             ini.C = c(-1, 0.1, 0.5, 1), 
             strip.num = 2000, control = list(), 
             fig.opt = TRUE, xlim = NULL,  
             ylim = NULL, unit=NULL, main=NULL){

  if(length(x)!=length(y)) 
      stop("'x' should have the same data length as 'y'!")
  Tem <- cbind(x, y)
  Tem <- na.omit(Tem)
  x   <- Tem[,1]
  y   <- Tem[,2]  

  poly1      <- owin(xrange=c(min(x)[1], max(x)[1]), 
                     yrange=c(min(y)[1],max(y)[1]))
  poly0      <- as.polygonal( owin(poly=list(x=x, y=y )) )
  data0      <- ppp(x, y, window= poly1) 

  area0      <- area.owin(poly0)
  perimeter0 <- perimeter(poly0)
  mat        <- pairdist(data0) 
    
  # To find the points associated with the egg tip and egg base
  qq <- cbind( which(mat == max(mat)[1], arr.ind=TRUE) )
  q  <- as.numeric( qq ) 
  x1 <- x[q[1]]
  y1 <- y[q[1]]
  x2 <- x[q[2]]
  y2 <- y[q[2]] 



  if( is.null(angle) ){
    xx      <- x-x2
    yy      <- y-y2
    x1      <- x1-x2
    y1      <- y1-y2
    x2      <- x2-x2
    y2      <- y2-y2
    z       <- sqrt( (y2-y1)^2 + (x2-x1)^2 )
    theta   <- acos( (x1-x2)/z ) %% (2*pi)
  }
  
  if( !is.null(angle) ){ 
    if( is.null(x0) | is.null(y0) )
        stop("When the angle is not null, x0 and y0 should be not null!")
    if( !is.numeric(angle) ) 
        stop("The argument angle should be a numerical value!")
    theta <- angle + pi
    xx    <- x-x0
    yy    <- y-y0
  }

  x.new   <- xx*cos(theta) + yy*sin(theta)
  y.new   <- yy*cos(theta) - xx*sin(theta)
  xmin    <- min(x.new)
  ymean   <- mean(y.new)
  x.new   <- x.new-xmin
  y.new   <- y.new-ymean 
  length0 <- max(x.new)[1]-min(x.new)[1]

  x.max   <- max(x.new)[1]
  x.min   <- min(x.new)[1]
  y.max   <- max(y.new)[1]
  y.min   <- min(y.new)[1]
  upper.win  <- owin( c(x.min, x.max), c(0, y.max) )
  lower.win  <- owin( c(x.min, x.max), c(y.min, 0) )
  total.poly <- as.polygonal(owin(poly=list(x=x.new, y=y.new)))
  poly1      <- intersect.owin(upper.win, total.poly)
  poly2      <- intersect.owin(lower.win, total.poly) 
  n          <- nrow(data.frame(total.poly))
  n1         <- nrow(data.frame(poly1))
  n2         <- nrow(data.frame(poly2)) 
  upper.area <- area.owin(poly1)
  lower.area <- area.owin(poly2)
  area0      <- area.owin( union.owin(poly1, poly2) )
  y.max      <- max( c(abs(y.max), abs(y.min)) )[1] 
  x.pred     <- seq( x.min, x.max, len=(strip.num+1) )
  part.upper.area <- c()
  part.lower.area <- c()
  widt            <- c()
  xmean           <- c()
  for(j in 2:(strip.num+1)){
    x.temp <- c(x.pred[j-1], x.pred[j])
    um.win <- owin( x.temp, c(0, y.max) ) 
    lm.win <- owin( x.temp, c(-y.max, 0) )
    m.area <- area.owin(um.win)
    u.conv <- setminus.owin(um.win, poly1)
    u.area <- m.area - area.owin(u.conv)
    l.conv <- setminus.owin(lm.win, poly2)
    l.area <- m.area - area.owin(l.conv)
    part.upper.area <- c(part.upper.area, u.area)
    part.lower.area <- c(part.lower.area, l.area)
    h1 <- NULL
    h2 <- NULL
    h1 <- try(intersect.owin(um.win, poly1), silent=TRUE)
    h2 <- try(intersect.owin(lm.win, poly2), silent=TRUE)
    g1 <- class(h1)
    g2 <- class(h2)
    if( g1=="owin" & g2=="owin" ){
        temp.u <- max(h1$y)[1]
        temp.l <- min(h2$y)[1]
        if(is.empty(h1)) temp.u <- 0
        if(is.empty(h2)) temp.l <- 0 
        widt   <- c(widt, temp.u-temp.l) 
        xmean  <- c(xmean, (x.pred[j-1] + x.pred[j])/2)   
    }
  }
  D        <- part.upper.area - part.lower.area
  xmean    <- xmean[!is.na(widt)]
  widt     <- widt[!is.na(widt)]
  width0   <- max(widt)[1]

  A0       <- length0
  B0       <- width0

  A1       <- A0/2
  A2       <- A0*3/4

  ind1     <- which.max(widt)
  C0       <- abs(xmean[ind1] - A0/2)
  temp2    <- abs(A2-xmean)
  ind2     <- which.min(temp2)
  D0       <- widt[ind2]

  ini.val  <- ini.C
  mat      <- matrix(NA, nrow=length(ini.val), ncol=2)

  obj.fun <- function(v){
      Par0  <- c(A0, B0, v, D0)
      z0    <- x.new - A0/2
      Ind1  <- which(y.new >= 0)
      Ind2  <- which(y.new <  0)
      xv1   <- z0[Ind1]
      yv1   <- y.new[Ind1]
      ypre1 <- NRGE(Par0, xv1)
      xv2   <- z0[Ind2]
      yv2   <- y.new[Ind2]
      ypre2 <- -NRGE(Par0, xv2)
      RSS0  <- sum( (yv1 - ypre1)^2 ) + sum( (yv2 - ypre2)^2 )
      return(RSS0)
  }
   
  for (i in 1:length(ini.val)) {
      res <- optim(ini.val[i], obj.fun, control=control)
      mat[i, ] <- c(res$par, res$val)
  }
 
  colnames(mat) <- c("w", "RSS")
  ind  <- which(mat[, 2] == min(mat[, 2])[1])[1]
  par  <- as.vector(mat[ind, 1])
  C0   <- par[1]

  Par0  <- c(A0, B0, C0, D0)
  z0    <- x.new - A0/2

  Ind1  <- which(y.new >= 0)
  Ind2  <- which(y.new <  0)
  xv1   <- z0[Ind1]
  yv1   <- y.new[Ind1]
  ypre1 <- NRGE(P=Par0, x=xv1)
  xv2   <- z0[Ind2]
  yv2   <- y.new[Ind2]
  ypre2 <- -NRGE(P=Par0, x=xv2)

  Ind3  <- sort(xv1, decreasing=TRUE, index.return=TRUE)$ix
  Ind4  <- sort(xv2, decreasing=FALSE, index.return=TRUE)$ix

  x.obs <- c(xv1[Ind3], xv2[Ind4])
  y.obs <- c(yv1[Ind3], yv2[Ind4])
  y.pre <- c(ypre1[Ind3], ypre2[Ind4])


  RSS0  <- sum( (y.obs-y.pre)^2 )
  RMSE0 <- sqrt( RSS0 / length(y.obs) )

  Resu  <- curveNRGE(P=c(0, 0, 0, Par0), 
      x=seq(-A0/2, A0/2, len=2000), fig.opt=F)

  x.theor <- Resu$x
  y.theor <- Resu$y

  if(!is.null(unit)){
    xlabel <- bquote( paste(italic("x"), " (", .(unit), ")", sep="") ) 
    ylabel <- bquote( paste(italic("y"), " (", .(unit), ")", sep="") )
  }
  if(is.null(unit)){
    xlabel <- bquote( italic("x") ) 
    ylabel <- bquote( italic("y") )
  }

  if(is.null(xlim)) xlim <- NULL
  if(is.null(ylim)) ylim <- NULL

  if(!is.null(xlim)) xlim <- xlim
  if(!is.null(ylim)) ylim <- ylim

  if(fig.opt == "TRUE"|fig.opt == "T"|fig.opt == "True"){ 
    dev.new()
    plot(z0, y.new, asp=1, col="grey50", lwd=2, 
         xlab=xlabel, ylab=ylabel, xlim=xlim, ylim=ylim,
         cex.lab=1.5, cex.axis=1.5, type="l")
    lines(x.theor, y.theor, col=2)
    title(main=main, cex.main=1.5, col.main=4, font.main=1)
    text(0, 0, bquote(paste("RMSE = ",  
      .(round(RMSE0,4)),sep="")), pos=NULL, cex=1.5, col=1)
  }

  list( theta=theta, x.obs=x.obs, y.obs=y.obs, y.pred=y.pre,
    par=c(A0,B0,C0,D0), scan.length=length0, scan.width=width0,   
    scan.area=area0, scan.perimeter=perimeter0, RSS=RSS0, 
    sample.size=length(y.new), RMSE=sqrt(RSS0/length(y.new)))

}

