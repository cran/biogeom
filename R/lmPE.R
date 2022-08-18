lmPE <- function(x, y, simpver = NULL, angle = NULL,              
             x0 = NULL, y0 = NULL, strip.num = 2000,  
             weights = NULL, fig.opt = TRUE, xlim = NULL,  
             ylim = NULL, unit = NULL, main = NULL){

  w <- weights
  if(!is.null(w) && !is.numeric(w)) 
     stop("'weights' must be a numeric vector")
  if(is.null(simpver)){
    if(!is.null(w) && length(w)!=4) 
      stop("'weights' should have four numeric elements")     
  }
  if(!is.null(simpver) && !is.numeric(simpver))
     stop("'simpver' must be chosen from 1, 2 and 3") 
  if(!is.null(simpver) && simpver==1){    
    if(!is.null(w) && length(w)!=3) 
      stop("'weights' should have three numeric elements")     
  }
  if(!is.null(simpver) && simpver==2){    
    if(!is.null(w) && length(w)!=2) 
      stop("'weights' should have two numeric elements")     
  }
  if(!is.null(simpver) && simpver==3){    
    if(!is.null(w) && length(w)!=2) 
      stop("'weights' should have two numeric elements")     
  }

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
  qq <- cbind( which(mat == max(mat)[1], arr.ind=TRUE)[1,] )
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
    theta <- angle - pi/2
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
  D      <- part.upper.area - part.lower.area
  xmean  <- xmean[!is.na(widt)]
  widt   <- widt[!is.na(widt)]
  width0 <- max(widt)[1]

  x.new  <- x.new-(max(x.new)[1]+min(x.new)[1])/2
  xv     <- x.new/(length0/2)
  yv     <- y.new/(length0/2)
  Index2 <- which(yv < 0)

  temp       <- 1-xv^2
  ind0       <- which(temp < 0)
  temp[ind0] <- 0
  z0         <- sqrt(temp)
  z1         <- xv*z0
  z2         <- xv^2*z0
  z3         <- xv^3*z0

  z0[Index2] <- -z0[Index2]   
  z1[Index2] <- xv[Index2]*z0[Index2]
  z2[Index2] <- xv[Index2]^2*z0[Index2]
  z3[Index2] <- xv[Index2]^3*z0[Index2]

  if(is.null(simpver)){
    res <- lm(yv~-1 + z0 + z1 + z2 + z3, weights=w)
  }
  if(!is.null(simpver)){
    if(simpver==1)
        res <- lm(yv~-1 + z0 + z1 + z2, weights=w)    
    if(simpver==2)
        res <- lm(yv~-1 + z0 + z1, weights=w) 
    if(simpver==3)
        res <- lm(yv~-1 + z0 + z2, weights=w) 
  }

  yv.pred         <- TSE(P=as.numeric(res$coefficients), x=xv, simpver=simpver)
  yv.pred[Index2] <- -yv.pred[Index2]
  y.pred          <- yv.pred*length0/2
  RSS1            <- sum((yv-yv.pred)^2)
  RSS2            <- sum((y.new-y.pred)^2)
  RMSE1           <- sqrt(sum((yv-yv.pred)^2)/length(yv))
  RMSE2           <- sqrt(sum((y.new-y.pred)^2)/length(y.new))

  if(!is.null(unit)){
    xlabel <- bquote( paste(italic(x), " (", .(unit), ")", sep="") ) 
    ylabel <- bquote( paste(italic(y), " (", .(unit), ")", sep="") )
  }
  if(is.null(unit)){
    xlabel <- bquote( italic(x) ) 
    ylabel <- bquote( italic(y) )
  }

  if(is.null(xlim)) xlim <- NULL
  if(is.null(ylim)) ylim <- NULL

  if(!is.null(xlim)) xlim <- xlim
  if(!is.null(ylim)) ylim <- ylim

  if(fig.opt == "TRUE"|fig.opt == "T"|fig.opt == "True"){ 

    dev.new()
    par(family="serif")
    par(mar=c(5,5,2,2))
    plot(xv, yv, asp=1, col="grey50", lwd=2, 
         xlab=bquote(italic(x)), ylab=bquote(italic(y)), xlim=xlim, ylim=ylim,
         cex.lab=1.5, cex.axis=1.5, type="l")
    lines(xv, yv.pred, col=2)
    title(main=main, cex.main=1.5, col.main=4, font.main=1)

    text(0, 0, bquote(paste("RMSE = ",  
      .(round(RMSE1,4)),sep="")), pos=NULL, cex=1.5, col=1)

    dev.new()
    par(family="serif")
    par(mar=c(5,5,2,2))
    plot(x.new, y.new, asp=1, col="grey50", lwd=2, 
         xlab=xlabel, ylab=ylabel, xlim=xlim, ylim=ylim,
         cex.lab=1.5, cex.axis=1.5, type="l")
    lines(x.new, y.pred, col=2)
    title(main=main, cex.main=1.5, col.main=4, font.main=1)
    text(0, 0, bquote(paste("RMSE = ",  
      .(round(RMSE2,4)),sep="")), pos=NULL, cex=1.5, col=1)
  }

  list( lm.tse=res, par=res$coefficients, theta=theta, 
        x.obs=x.new, y.obs=y.new, y.pred=y.pred, x.stand.obs=xv, 
        y.stand.obs=yv, y.stand.pred=yv.pred, scan.length=length0, 
        scan.width=width0, scan.area=area0, scan.perimeter=perimeter0, 
        RSS=RSS2, sample.size=length(y.new), RMSE=RMSE2 )
}


