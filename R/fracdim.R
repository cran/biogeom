fracdim <- function(x, y,  frac.fig = TRUE,
             denomi.range = seq(8, 30, by = 1), 
             ratiox = 0.02, ratioy = 0.08, main = NULL){

    if(length(x)!=length(y)) 
        stop("'x' should have the same data length as 'y'!")
    Tem <- cbind(x, y)
    Tem <- na.omit(Tem)
    x   <- Tem[,1]
    y   <- Tem[,2]  

    poly1 <- owin(xrange=c(min(x)[1], max(x)[1]), 
                  yrange=c(min(y)[1],max(y)[1]))
    data0 <- ppp(x, y, window=poly1)
    xj    <- x
    yj    <- y 
    # Calculating a fractal dimension
    x.dis <- range(x)[2]-range(x)[1]
    y.dis <- range(y)[2]-range(y)[1]
    L     <- max(x.dis, y.dis)[1]
    xt    <- xj - min(xj)[1] 
    yt    <- yj - min(yj)[1]
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
    if(C2-C1 < 4.0) ylim0 <- c(C2-4.0, C2)


    xloc   <- xlim0[1] + (xlim0[2]-xlim0[1])*ratiox
    yloc   <- ylim0[2] - (ylim0[2]-ylim0[1])*ratioy
    if( frac.fig =="TRUE" | frac.fig == "T" | frac.fig == "True"){
      dev.new()
      plot( log.x, log.y, cex.lab=1.5, cex.axis=1.5, cex=1.5, xaxs="i", yaxs="i",
            xlim=xlim0, ylim=ylim0, las=1,
            xlab=expression(paste("ln ", delta, {""}^{-"1"}, sep="")), 
            ylab=expression(paste("ln ", italic("N"), sep="")) )
      title(main=main, cex.main=1.5, col.main=4, font.main=1)
      abline( resu, col=2 )   
      text(xloc, yloc-0*(ylim0[2]-ylim0[1])/10, bquote( paste(italic("y = "),  
        .(round(a, 3)), " + ", .(round(b, 3)), italic("x"), seq="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-1*(ylim0[2]-ylim0[1])/10, bquote(paste("CI: ", 
        .(round(confint( resu )[2,1], 3)), ", ", .(round(confint( resu )[2,2], 3)),  
        sep="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-2*(ylim0[2]-ylim0[1])/10, 
        bquote( paste(italic(r), {""}^{"2"}, " = ", 
        .(round(summary(resu)$r.squared, 3)), sep="") ), pos=4, cex=1.5, col=1)   
      text(xloc, yloc-3*(ylim0[2]-ylim0[1])/10, 
        bquote( paste(italic(n), , " = ", .(length(log.x)),
        sep="") ), pos=4, cex=1.5, col=1) 
    }
    list( a=a, sd.a=sd.a, lci.a=lci.a, uci.a=uci.a, b=b, sd.b=sd.b, 
          lci.b=lci.b, uci.b=uci.b, r.sq=r.sq, delta=delta, N=N  ) 
  }