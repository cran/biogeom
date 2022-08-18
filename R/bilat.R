bilat <- function(x, y, strip.num = 200, peri.np = NULL, n.loop = 60,
           auto.search = TRUE, animation.fig = TRUE, time.interval = 0.001,  
           unit = "cm", main = NULL, diff.fig = TRUE,
           angle = NULL, ratiox = 0.02, ratioy = 0.08, 
           fd.opt = TRUE, frac.fig = TRUE,
           denomi.range = seq(8, 30, by = 1)){

  if( !is.null(angle) & !is.numeric(angle) ) 
      stop("'angle' should be a numerical value if it is not null.")

  if(length(x)!=length(y)) 
      stop("'x' should have the same data length as 'y'!")

  Tem <- cbind(x, y)
  Tem <- na.omit(Tem)
  x   <- Tem[,1]
  y   <- Tem[,2]  

  if(strip.num / floor(strip.num) != 1 | strip.num <= 0){
    stop("'strip.num' must be a positive integer")
  }

  if(is.null(peri.np)){
    poly0 <- as.polygonal( owin(poly=list(x=x, y=y)) )
    perimeter0 <- perimeter(poly0)
  }

  if(!is.null(peri.np)){
    if(!is.numeric(peri.np)){
      stop("'peri.np' must be a poitive integer!")
    }
    if(is.numeric(peri.np) & peri.np <= 0){
      stop("'peri.np' must be a poitive integer!")
    }
    my.num <- length(x)
    if( length(x) > peri.np ) my.num <- peri.np 
  
    my.peri <- c()  
    my.ind  <- sort( sample(1:length(x), n.loop, replace=FALSE) )
    for(i in 1:n.loop){
      if(my.ind[i] == 1){
         x1  <- x[1:length(x)]
         y1  <- y[1:length(x)]
         ind <- floor(as.numeric(quantile(1:length(x1), 
                      probs=seq(0, 1, len=my.num))))
         x2  <- x1[ind]
         y2  <- y1[ind]             
      }
      if(my.ind[i] > 1){
        inde <- c(my.ind[i]:length(x), 1:(my.ind[i]-1))
        x1   <- x[inde]
        y1   <- y[inde]  
        ind  <- floor(as.numeric(quantile(1:length(x1), 
                      probs=seq(0, 1, len=my.num))))
        x2   <- x1[ind]
        y2   <- y1[ind]    
      }
      my.poly <- as.polygonal( owin(poly=list(x=x2, y=y2)) )
      my.peri <- c(my.peri, perimeter(my.poly))  
    }
    perimeter0 <- mean(my.peri)   
  }

  poly1      <- owin(xrange=c(min(x)[1], max(x)[1]), 
                     yrange=c(min(y)[1],max(y)[1]))
  data0      <- ppp(x, y, window=poly1)
  x          <- data0$x
  y          <- data0$y
  n          <- length(x)
  mat        <- pairdist(data0) 
  xlabel <- bquote(paste(italic(x)," (", .(unit),")",sep="")) 
  ylabel <- bquote(paste(italic(y)," (", .(unit),")",sep=""))
  if(auto.search != "TRUE" & auto.search != "T"){
    dev.new()
    plot( x, y, asp=1, type="l", xlab=xlabel, ylab=ylabel, 
          cex.lab=1.5, cex.axis=1.5, lwd=3, col="grey50" )
    resu  <- locator(n = 2, type="o", col=2, cex=1.5)
    #### (x1, y1) represents the location of the leaf apex
    #### (x2, y2) represents the location of the leaf base
    x1    <- resu$x[1]
    x2    <- resu$x[2]
    y1    <- resu$y[1]
    y2    <- resu$y[2]
    if(x1 < x2 | y1 < y2){
stop("\n Leaf apex must be on the top right corner of leaf base! 
\n The angle from leaf apex to base should range 0 to pi/2!
\n The leaf apex should be chosen earlier than the leaf base!")
    }
  #### To find a point on the edge that has the shortest 
  ####   distance from a chosen point
  #### Two points on the edge will be found 
  ####   corresponding to the leaf apex and base
    dis.fun1 <- function(z){
      sqrt((z[1]-x1)^2 + (z[2]-y1)^2)    }
    pp       <- cbind(x,y)
    dis.val1 <- apply(pp, 1L, dis.fun1)
    ind1     <- which(dis.val1==min(dis.val1)[1])[1]
    x1       <- pp[ind1,1]
    y1       <- pp[ind1,2]
    dis.fun2 <- function(z){
      sqrt((z[1]-x2)^2 + (z[2]-y2)^2)
    }
    dis.val2 <- apply(pp, 1L, dis.fun2)
    ind2     <- which(dis.val2==min(dis.val2)[1])[1]
    x2       <- pp[ind2,1]
    y2       <- pp[ind2,2]

  ####################################################################
  }
  if(auto.search == "TRUE" | auto.search == "T"){
    # To find the points corresponding 
    #   to the leaf apex and leaf base
    qq <- cbind( which(mat == max(mat)[1], arr.ind=TRUE)[1,] )
    q  <- as.numeric( qq ) 
    x1 <- x[q[1]]
    y1 <- y[q[1]]
    x2 <- x[q[2]]
    y2 <- y[q[2]] 
  } 
  length0 <- sqrt((x1-x2)^2+(y1-y2)^2)
  length0 <- length0[[1]]
  x       <- x-x2
  y       <- y-y2
  x1      <- x1-x2
  y1      <- y1-y2
  x2      <- x2-x2
  y2      <- y2-y2
  z       <- sqrt( (y2-y1)^2 + (x2-x1)^2 )
  theta   <- acos( (x1-x2)/z ) %% (2*pi)
  x.new   <- x*cos(theta) + y*sin(theta)
  y.new   <- y*cos(theta) - x*sin(theta)
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
        xmean  <- c( xmean, (x.pred[j-1] + x.pred[j])/2 )     
    }
  }
  D        <- part.upper.area - part.lower.area
  xmean    <- xmean[!is.na(widt)] # Former
  widt     <- widt[!is.na(widt)]  # Latter
  max.ind  <- which.max(widt)[1]
  width0   <- widt[max.ind]
  x.width  <- xmean[max.ind]
  x.1e     <- length0*1/8
  x.2e     <- length0*2/8
  x.4e     <- length0*4/8
  x.6e     <- length0*6/8
  x.7e     <- length0*7/8
  temp1    <- abs(x.1e-xmean)
  IND1     <- which.min(temp1)
  width.1e <- widt[IND1]
  temp2    <- abs(x.2e-xmean)
  IND2     <- which.min(temp2)
  width.2e <- widt[IND2]
  temp3    <- abs(x.4e-xmean)
  IND3     <- which.min(temp3)
  width.4e <- widt[IND3]
  temp4    <- abs(x.6e-xmean)
  IND4     <- which.min(temp4)
  width.6e <- widt[IND4]
  temp5    <- abs(x.7e-xmean)
  IND5     <- which.min(temp5)
  width.7e <- widt[IND5]
  bi.test  <- wilcox.test(part.upper.area, 
                part.lower.area, paired=TRUE)
  abs.D    <- abs(D)
  SI       <- sum(abs.D/(part.upper.area+
                part.lower.area))/strip.num
  AR       <- upper.area / lower.area
  Xlab     <- bquote(paste("Difference (",
                .(unit), ''^'2'*')', sep="")) 

  if(is.null(angle))  phi <- theta
  if(!is.null(angle)) phi <- angle

  if(animation.fig =="TRUE"|animation.fig==
          "T"|animation.fig == "True"){
    if(phi %% (2*pi) != 0){
      xx.new  <- x.new*cos(phi) - y.new*sin(phi)
      yy.new  <- y.new*cos(phi) + x.new*sin(phi)   
    }
    else{
      xx.new  <- x.new
      yy.new  <- y.new
    }
    dev.new()
    plot( xx.new, yy.new, type="n", 
          xlab=xlabel, ylab=ylabel, cex=0.1, pch=16,
          cex.lab=1.5, cex.axis=1.5, asp=1)
    title(main=main, cex.main=1.5, col.main=4, font.main=1)
    for(i in 1:length(xx.new)){
      Sys.sleep(time.interval)
      points(xx.new[i], yy.new[i], cex=0.1, col="grey70", pch=16)
    }
    polygon(xx.new, yy.new, col = "grey90", border = "grey70")
    abline(v=0, col=4, lwd=1, lty=2)
    abline(h=0, col=4, lwd=1, lty=2)
    if(phi %% (pi/2) != 0 ){
      abline(0, tan(phi), col=2, lwd=1, lty=2)
    }
    lines(xx.new, yy.new, lwd=1, col="grey70")
  }

  if(diff.fig =="TRUE" | diff.fig == "T" | diff.fig == "True"){
      dev.new()     
      colors <- ifelse(D >= 0, "blue", "red")
      plot( D, type="h", ylab=Xlab, xlab="Strip number", 
            cex.lab=1.5, cex.axis=1.5, lwd=1, col=colors )
      abline(h=0, col="black", lwd=1)
  }

  if(fd.opt=="TRUE"|fd.opt=="T"|fd.opt=="True"){
    poly0 <- as.polygonal( owin(poly=list(x=x, y=y)) )
    xj    <- x
    yj    <- y 
    # Calculating a fractal dimension of the polygon's perimeter
    x.dis <- range(x)[2]-range(x)[1]
    y.dis <- range(y)[2]-range(y)[1]
    L     <- max(x.dis, y.dis)[1]
    xt    <- xj - min(xj)[1] 
    yt    <- yj - min(yj)[1]
    #xt    <- xt + (z-max(xt)[1])/2
    #yt    <- yt + (z-max(yt)[1])/2
    Side     <- L/denomi.range
    k       <- 0
    delta   <- c()
    N       <- c()
    for(i in 1:length(Side)){
      k      <- k + 1
      ndiv   <- L/Side[i]
      W      <- owin(c(0, L), c(0, L))
      Index  <- inside.owin(xt, yt, W)
      xt2    <- xt[Index]
      yt2    <- yt[Index]     
      Z      <- ppp(xt2, yt2, window=W)
      Number <- array( quadratcount(Z, ndiv, ndiv) ) 
      N      <- c(N, length(Number[Number > 0]))
      delta  <- c(delta, Side[i])
    }
    log.x  <- -log(delta)
    log.y  <- log(N)
    resu   <- lm(log.y~log.x)    
    a      <- summary(resu)$coefficients[[1]]
    b      <- summary(resu)$coefficients[[2]]
    sd.a   <- summary(resu)$coefficients[[3]]
    sd.b   <- summary(resu)$coefficients[[4]]
    lci.a  <- confint( resu )[1,1]
    uci.a  <- confint( resu )[1,2]
    lci.b  <- confint( resu )[2,1]
    uci.b  <- confint( resu )[2,2]
    r.sq   <- summary(resu)$r.squared
    x0     <- log.x
    y0     <- log.y
    B  <- min(x0)[1]
    u  <- NA
    for(s in seq(-100, 100, by=0.5)){
      if(B[1] < s + 0.5 & B[1] >= s) u <- s
    }
    B1 <- u
    B  <- max(x0)[1]
    u  <- NA
    for(s in seq(-100, 100, by=0.5)){
      if(B[1] < s + 0.5 & B[1] >= s) u <- s + 0.5
    }
    B2 <- u
    C  <- min(y0)[1]
    u  <- NA
    for(s in seq(-100, 100, by=1)){
      if(C[1] < s + 1 & C[1] >= s) u <- s
    }
    C1 <- u
    C  <- max(y0)[1]
    u  <- NA
    for(s in seq(-100, 100, by=1)){
      if(C[1] < s + 1 & C[1] >= s) u <- s + 1
    }
    C2 <- u
    xlim0    <- c(B1, B2)
    ylim0    <- c(C1, C2) 
    if(B2-B1 < 2.0) xlim0 <- c(B2-2.0, B2)
    # if(C2-C1 < 4.0) ylim0 <- c(C2-4.0, C2)     
    if(C2-C1 < 2.5) ylim0 <- c(C2-2.5, C2)
    
    xloc   <- xlim0[1] + (xlim0[2]-xlim0[1])*ratiox
    yloc   <- ylim0[2] - (ylim0[2]-ylim0[1])*ratioy
    if( frac.fig =="TRUE" | frac.fig == "T" | frac.fig == "True"){
      dev.new()
      plot( log.x, log.y, cex.lab=1.5, cex.axis=1.5, 
            cex=1.5, xaxs="i", yaxs="i",
            xlim=xlim0, ylim=ylim0, las=1,
            xlab=expression(paste("ln ", 
            italic("\u03B4"), {""}^{-"1"}, sep="")), 
            ylab=expression(paste("ln ", italic("N"), sep="")) )
      title(main=main, cex.main=1.5, col.main=4, font.main=1)
      abline( resu, col=2 )   
      text(xloc, yloc-0*(ylim0[2]-ylim0[1])/10, 
           bquote( paste(italic("y = "),  
           .(round(a, 3)), " + ", .(round(b, 3)), 
           italic("x"), seq="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-1*(ylim0[2]-ylim0[1])/10, 
           bquote(paste("CI: ", 
           .(round(confint( resu )[2,1], 3)), ", ", 
           .(round(confint( resu )[2,2], 3)),  
           sep="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-2*(ylim0[2]-ylim0[1])/10, 
        bquote( paste(italic(r), {""}^{"2"}, " = ", 
        .(round(summary(resu)$r.squared, 3)), sep="") ), 
        pos=4, cex=1.5, col=1)   
      text(xloc, yloc-3*(ylim0[2]-ylim0[1])/10, 
        bquote( paste(italic(n), , " = ", .(length(log.x)),
        sep="") ), pos=4, cex=1.5, col=1) 
    }
    te <- list( x=xx.new, y=yy.new, phi=phi, n1=n1, n2=n2, n=n, 
          total.poly=total.poly, upper.poly=poly1, 
          lower.poly=poly2, part.upper.area=part.upper.area,  
          part.lower.area=part.lower.area, D=D,   
          upper.area=upper.area, lower.area=lower.area, 
          SI=SI, AR=AR, scan.length=length0, 
          scan.width=width0, scan.area=area0, 
          scan.perimeter=perimeter0, x.width=x.width, width.1e=width.1e,
          width.2e=width.2e, width.4e=width.4e, width.6e=width.6e,
          width.7e=width.7e, bi.test=bi.test,
          a=a, sd.a=sd.a, lci.a=lci.a, uci.a=uci.a, b=b, sd.b=sd.b, 
          lci.b=lci.b, uci.b=uci.b, r.sq=r.sq, delta=delta, N=N  ) 
  }
  if(fd.opt != "TRUE" & fd.opt != "T" & fd.opt != "True"){
    te <- list( x=xx.new, y=yy.new, phi=phi, n1=n1, n2=n2, n=n, 
          total.poly=total.poly, upper.poly=poly1, 
          lower.poly=poly2, part.upper.area=part.upper.area,  
          part.lower.area=part.lower.area, D=D,   
          upper.area=upper.area, lower.area=lower.area, 
          SI=SI, AR=AR, scan.length=length0, 
          scan.width=width0, scan.area=area0, 
          scan.perimeter=perimeter0, x.width=x.width, width.1e=width.1e,
          width.2e=width.2e, width.4e=width.4e, width.6e=width.6e,
          width.7e=width.7e, bi.test=bi.test )
  }
  
  return( te )
}

