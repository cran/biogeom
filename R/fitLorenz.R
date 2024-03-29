fitLorenz <- function(expr, z, ini.val, simpver = 4, 
               control = list(), par.list = FALSE, 
               fig.opt = FALSE, np = 2000, 
               xlab=NULL, ylab=NULL, main = NULL, subdivisions = 100L,
               rel.tol = .Machine$double.eps^0.25, 
               abs.tol = rel.tol, stop.on.error = TRUE, 
               keep.xy = FALSE, aux = NULL, par.limit = TRUE){


  if(!(par.limit %in% c("F", "False", "FALSE", "T", "True", "TRUE")))
    stop("'par.limit' only can be set to be FALSE or TRUE!")

  # The code for fitting the original Lorenz functions
  if(!identical(expr, MPerformanceE)){

    if(length(z) < 2) 
        stop("'z' should include 2 elements at least!")

    z   <- na.omit(z)
    z1  <- sort( z )
    x1  <- 1:length(z1)/length(z1)
    y1  <- cumsum(z1)/sum(z1)

    # Rotates the data and then shifts them to the right by a distance of sqrt(2) 
    theta <- -3/4*pi
    x     <- x1 * cos(theta) + y1 * sin(theta)
    y     <- y1 * cos(theta) - x1 * sin(theta)
    x     <- x+sqrt(2)

    ini.val <- as.list(ini.val)
    p       <- length(ini.val)
    s       <- 1
    for (i in 1:p) {
        s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat     <- matrix(NA, nrow = s, ncol = (p + 1))

    if(par.limit == "T" | par.limit == "TRUE" | par.limit == "True"){
      if(identical(expr, SHE)){
        obj.fun <- function(P){
          if(P[1] < 0 | P[1] >= 1 | P[2] < 0 | P[2] > 1 | P[3] < 0 | P[3] > 1 | P[4] < 1)
              RSS <- Inf
          if(P[1] >= 0 & P[1] < 1 & P[2] >= 0 & P[2] <= 1 & P[3] >= 0 & P[3] <= 1 & P[4] >= 1){
            y1.theo <- expr(P, x1)
            theta   <- -3/4*pi
            x       <- x1 * cos(theta) + y1 * sin(theta)
            y.theo  <- y1.theo * cos(theta) - x1 * sin(theta)
            x       <- x+sqrt(2)
            y       <- y1 * cos(theta) - x1 * sin(theta)
            RSS    <- sum((y.theo-y)^2) 
          }
          return(RSS) 
        }
      }

      if(identical(expr, SCSE)){
        obj.fun <- function(P){
          if(P[1] < 0 | P[2] <= 0 | P[2] > 1 | P[3] < 1)
              RSS <- Inf
          if(P[1] >= 0 & P[2] > 0 & P[2] <= 1 & P[3] >= 1){
              y1.theo <- expr(P, x1)
              theta   <- -3/4*pi
              x       <- x1 * cos(theta) + y1 * sin(theta)
              y.theo  <- y1.theo * cos(theta) - x1 * sin(theta)
              x       <- x+sqrt(2)
              y       <- y1 * cos(theta) - x1 * sin(theta) 
              RSS    <- sum((y.theo-y)^2) 
          }
          return(RSS) 
        }
      }

      if(identical(expr, SarabiaE)){
        obj.fun <- function(P){
          if(P[1] < 0 | P[2]*P[4] + P[1] > 1 | P[2]*P[4] < 0 | P[3] < 0 | P[2] + 1 <= 0)
              RSS <- Inf
          if(P[1] >= 0 & P[2]*P[4] + P[1] <= 1 & P[2]*P[4] >= 0 & P[3] >= 0 & P[2] + 1 > 0){
            y1.theo <- expr(P, x1)
            theta   <- -3/4*pi
            x       <- x1 * cos(theta) + y1 * sin(theta)
            y.theo  <- y1.theo * cos(theta) - x1 * sin(theta)
            x       <- x+sqrt(2)
            y       <- y1 * cos(theta) - x1 * sin(theta) 
            RSS    <- sum((y.theo-y)^2) 
          }
          return(RSS) 
        }
      }
    }

    if(par.limit == "F" | par.limit == "FALSE" | par.limit == "False"){
      obj.fun <- function(P){
          y1.theo <- expr(P, x1)
          theta   <- -3/4*pi
          x       <- x1 * cos(theta) + y1 * sin(theta)
          y.theo  <- y1.theo * cos(theta) - x1 * sin(theta)
          x       <- x+sqrt(2)
          y       <- y1 * cos(theta) - x1 * sin(theta) 
          RSS     <- sum((y.theo-y)^2) 
      }
    }

    for (i in 1:nrow(ini.val)) {
        res <- optim(ini.val[i, ], obj.fun, control = control)
        mat[i, ] <- c(res$par, res$val)
    }
    Names <- rep(NA, len = p)
    for (k in 1:p) {
        Names[k] <- paste("P[", k, "]", sep = "")
    }
    colnames(mat) <- c(Names, "RSS")

    ind    <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
    par    <- as.vector(mat[ind, 1:p])

    x3      <- seq(0, 1, len=np)
    y3      <- expr(par, x3)
    yy2     <- expr(par, x1)
    theta   <- -3/4*pi
    x2      <- x3 * cos(theta) + y3 * sin(theta)   
    y2      <- y3 * cos(theta) - x3 * sin(theta)  
    x2      <- x2 + sqrt(2)                        
    y.theo  <- yy2 * cos(theta) - x1 * sin(theta)  

    RSS    <- sum((y.theo - y)^2)
    r.sq   <- 1-sum((y.theo - y)^2)/sum((y-mean(y))^2)

    goal.fun <- function(x){
      (x - expr(par, x))*2
    }  

    if(identical(expr, SHE)){
      GC <- 1-2*(1-par[1])/(par[4]+1)
    }

    if(identical(expr, SarabiaE)){
      GC <- par[1]*(1-2/(2+par[3])) + par[2]*(1-2/(2+par[4]))
    }

    if(!identical(expr, SHE) & !identical(expr, SarabiaE)){
      GC <- integrate(goal.fun, lower=0, upper=1, 
            subdivisions = subdivisions, 
            rel.tol = rel.tol, abs.tol = abs.tol, 
            stop.on.error = stop.on.error, keep.xy = keep.xy, 
            aux = aux)$value
    }

    if(is.null(xlab))
       xlab <- "Cumulative proportion of number"
    if(is.null(ylab))
       ylab <- "Cumulative proportion of size"

    if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True"){
      ratiox <- 0.02
      ratioy <- 0.08
      xlim1  <- c(0, 1)
      ylim1  <- c(0, 1) 
      xloc   <- xlim1[1] + (xlim1[2]-xlim1[1])*ratiox
      yloc   <- ylim1[2] - (ylim1[2]-ylim1[1])*ratioy
      xlim2  <- c(0, 1.5)
      y2.candi <- as.numeric( na.omit(y2) )
      ymax <- max(c(y, y2.candi))[1]   
      A1   <- 10^seq(-6, 1, by=1)
      B1   <- A1*0.08
      B2   <- A1*0.12
      B3   <- A1*0.16
      B4   <- A1*0.20
      B5   <- A1*0.25
      B6   <- A1*0.40
      B7   <- A1*0.50
      B8   <- A1*0.70
      B9   <- A1*0.80   

      for(i in 1:length(A1)){
        if(ymax < B1[i])                { ylim2 <- c(0, B1[i]); break }
        if(ymax >= B1[i] & ymax < B2[i]){ ylim2 <- c(0, B2[i]); break }
        if(ymax >= B2[i] & ymax < B3[i]){ ylim2 <- c(0, B3[i]); break }
        if(ymax >= B3[i] & ymax < B4[i]){ ylim2 <- c(0, B4[i]); break }
        if(ymax >= B4[i] & ymax < B5[i]){ ylim2 <- c(0, B5[i]); break }
        if(ymax >= B5[i] & ymax < B6[i]){ ylim2 <- c(0, B6[i]); break }
        if(ymax >= B6[i] & ymax < B7[i]){ ylim2 <- c(0, B7[i]); break }
        if(ymax >= B7[i] & ymax < B8[i]){ ylim2 <- c(0, B8[i]); break }
        if(ymax >= B8[i] & ymax < B9[i]){ ylim2 <- c(0, B9[i]); break }
        if(ymax >= B9[i] & ymax < A1[i]){ ylim2 <- c(0, A1[i]*1.2); break }
      }  

      dev.new()
      layout( matrix(1:2, 1, 2, byrow=TRUE) )

      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )
      plot( x1, y1, main=main, cex.lab=1.5, cex.axis=1.5, 
            las=1, xaxs="i", yaxs="i", 
            xlab=xlab, ylab=ylab, 
            xlim=xlim1, ylim=ylim1, type="n" )
      polygon(x3, y3, col="grey70")
      points(x1, y1, cex=1.5, pch=16)
      lines(x3, y3, col=2) 
      lines(c(0, 1), c(0, 1), col=1)

      if(identical(expr, SHE)){
        if(par[1] >= 1e-9)
          text(xloc, yloc, bquote(paste(hat(delta), " = ", 
            .(sprintf("%.4f", round(par[1], 4))), sep="")), pos=4, cex=1.5, col=1)
        if(par[1] < 1e-9)
          text(xloc, yloc, bquote(paste(hat(delta), " = ", 
            .(sprintf("%.0f", round(par[1], 0))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-1.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(rho), " = ", 
           .(sprintf("%.4f", round(par[2], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-2.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(omega), " = ", 
          .(sprintf("%.4f", round(par[3], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(P)), " = ", 
             .(sprintf("%.4f", round(par[4], 4))), sep="")), pos=4, cex=1.5, col=1) 
        text(xloc, yloc-4.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
             .(sprintf("%.4f", round(r.sq, 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-5.0*(ylim1[2]-ylim1[1])/10, 
             bquote( paste(italic(n), " = ", 
             .(length(x)), sep="") ), pos=4, cex=1.5, col=1)
      }

      if(identical(expr, SCSE)){
        if(par[1] >= 1e-9)
          text(xloc, yloc, bquote(paste(hat(gamma), " = ", 
            .(sprintf("%.4f", round(par[1], 4))), sep="")), pos=4, cex=1.5, col=1)
        if(par[1] < 1e-9)
          text(xloc, yloc, bquote(paste(hat(gamma), " = ", 
            .(sprintf("%.0f", round(par[1], 0))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-1.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(alpha), " = ", 
           .(sprintf("%.4f", round(par[2], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-2.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(beta), " = ", 
          .(sprintf("%.4f", round(par[3], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
             .(sprintf("%.4f", round(r.sq, 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-4.0*(ylim1[2]-ylim1[1])/10, 
             bquote( paste(italic(n), " = ", 
             .(length(x)), sep="") ), pos=4, cex=1.5, col=1)
      }

      if(identical(expr, SarabiaE)){
        text(xloc, yloc-0.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(lambda), " = ", 
           .(sprintf("%.4f", round(par[1], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-1.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(eta), " = ", 
          .(sprintf("%.4f", round(par[2], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-2.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(a)), {""}[1], " = ", 
             .(sprintf("%.4f", round(par[3], 4))), sep="")), pos=4, cex=1.5, col=1) 
        if(abs(par[4]) >= 1e-9)
          text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(a)), {""}[2], " = ", 
            .(sprintf("%.4f", round(par[4], 4))), sep="")), pos=4, cex=1.5, col=1)
        if(abs(par[4]) < 1e-9)
          text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(a)), {""}[2], " = ",
            .(sprintf("%.0f", round(par[4], 0))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-4.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
             .(sprintf("%.4f", round(r.sq, 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-5.0*(ylim1[2]-ylim1[1])/10, 
             bquote( paste(italic(n), " = ", 
             .(length(x)), sep="") ), pos=4, cex=1.5, col=1)
      }

      if(!identical(expr, SHE) & !identical(expr, SCSE) & !identical(expr, SarabiaE)){
        text(xloc, yloc-0.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
             .(sprintf("%.4f", round(r.sq, 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-1.0*(ylim1[2]-ylim1[1])/10, 
             bquote( paste(italic(n), " = ", 
             .(length(x)), sep="") ), pos=4, cex=1.5, col=1)
      }
   
      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )

      plot(x, y, main=main, cex.lab=1.5, cex.axis=1.5,   
           las=1, xaxs="i", yaxs="i", xlab=expression(italic(x)), 
           ylab=expression(italic(y)), xlim=xlim2, ylim=ylim2, type="n")
      points(x, y, cex=1.5, pch=16)
      lines(x2, y2, col=2, lty=1, lwd=1) 

    }

    para.tab <- data.frame( Parameter = c(Names, "r.sq", "RSS", 
                    "sample.size", "GC"), Estimate = c(par, r.sq, RSS, 
                    length(x), GC) )
    if(par.list == "T" | par.list == "TRUE" | par.list == "True"){
        print(para.tab)
        cat("\n")
    }
    return(list( x1=x1, y1=y1, x=x, y=y, par=par, r.sq=r.sq, 
      RSS=RSS, sample.size=length(x), GC=GC) )

  }




  # The code for fitting the rotated and right-shifted Lorenz functions, i.e., MPerformanceE
  if(identical(expr, MPerformanceE)){

    if(!identical(expr, MPerformanceE)){
      stop("'expr' should be 'MPerformanceE'")
    }

    if( !(simpver %in% c(4, 5)) )
    stop("'simpver' should be 4 or 5 when using 'MPerformanceE' to fit the rotated and right-shifed Lorenz curve!")  

    if(length(z) < 2) 
        stop("'z' should include 2 elements at least!")

    z   <- na.omit(z)
    z1  <- sort( z )
    x1  <- 1:length(z1)/length(z1)
    y1  <- cumsum(z1)/sum(z1)  

    # Rotates the data and then shifts them to the right by a distance of sqrt(2) 
    theta <- -3/4*pi
    x     <- x1 * cos(theta) + y1 * sin(theta)
    y     <- y1 * cos(theta) - x1 * sin(theta)
    x     <- x+sqrt(2)

    ini.val <- as.list(ini.val)
    p       <- length(ini.val)
    if(simpver==4 & p!=5)
        stop("There should be five parameter for version 4!")
    if(simpver==5 & p!=3)
        stop("There should be three parameter for version 5!")
    s       <- 1
    for (i in 1:p) {
        s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat     <- matrix(NA, nrow = s, ncol = (p + 1))
    obj.fun <- function(P){
      if(P[1] <= 0 | P[2] < 0 | P[3] < 0)
          RSS <- Inf
      if(P[1] > 0 | P[2] >= 0 | P[3] >= 0){
          y.theo <- expr(P=P, x=x, simpver=simpver)
          RSS    <- sum((y.theo-y)^2) 
      }
      return(RSS) 
    }
    for (i in 1:nrow(ini.val)) {
        res <- optim(ini.val[i, ], obj.fun, control = control)
        mat[i, ] <- c(res$par, res$val)
    }
    Names <- rep(NA, len = p)
    for (k in 1:p) {
        Names[k] <- paste("P[", k, "]", sep = "")
    }
    colnames(mat) <- c(Names, "RSS")

    ind    <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
    par    <- as.vector(mat[ind, 1:p])
    x2     <- seq(0, sqrt(2), len=np)
    y2     <- expr(P=par, x=x2, simpver=simpver)
    y.theo <- expr(P=par, x=x, simpver=simpver)
    RSS    <- sum((y.theo - y)^2)
    r.sq   <- 1-sum((y.theo - y)^2)/sum((y-mean(y))^2)
    goal.fun <- function(x){
      expr(P=par, x=x, simpver=simpver)
    }  
    GC <- integrate(goal.fun, lower=0, upper=sqrt(2), 
            subdivisions = subdivisions, 
            rel.tol = rel.tol, abs.tol = abs.tol, 
            stop.on.error = stop.on.error, keep.xy = keep.xy, 
            aux = aux)$value/0.5

    theta <- 3/4*pi
    xt    <- x2 - sqrt(2)
    x3    <- xt * cos(theta) + y2 * sin(theta)
    y3    <- y2 * cos(theta) - xt * sin(theta)


    if(is.null(xlab))
       xlab <- "Cumulative proportion of number"
    if(is.null(ylab))
       ylab <- "Cumulative proportion of size"

    if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True"){
      ratiox <- 0.02
      ratioy <- 0.08
      xlim1  <- c(0, 1)
      ylim1  <- c(0, 1) 
      xloc   <- xlim1[1] + (xlim1[2]-xlim1[1])*ratiox
      yloc   <- ylim1[2] - (ylim1[2]-ylim1[1])*ratioy
      xlim2  <- c(0, 1.5)

      ymax <- max(c(y, y2))[1]    
      A1   <- 10^seq(-6, 1, by=1)
      B1   <- A1*0.08
      B2   <- A1*0.12
      B3   <- A1*0.16
      B4   <- A1*0.20
      B5   <- A1*0.25
      B6   <- A1*0.40
      B7   <- A1*0.50
      B8   <- A1*0.70
      B9   <- A1*0.80   

      for(i in 1:length(A1)){
        if(ymax < B1[i])                { ylim2 <- c(0, B1[i]); break }
        if(ymax >= B1[i] & ymax < B2[i]){ ylim2 <- c(0, B2[i]); break }
        if(ymax >= B2[i] & ymax < B3[i]){ ylim2 <- c(0, B3[i]); break }
        if(ymax >= B3[i] & ymax < B4[i]){ ylim2 <- c(0, B4[i]); break }
        if(ymax >= B4[i] & ymax < B5[i]){ ylim2 <- c(0, B5[i]); break }
        if(ymax >= B5[i] & ymax < B6[i]){ ylim2 <- c(0, B6[i]); break }
        if(ymax >= B6[i] & ymax < B7[i]){ ylim2 <- c(0, B7[i]); break }
        if(ymax >= B7[i] & ymax < B8[i]){ ylim2 <- c(0, B8[i]); break }
        if(ymax >= B8[i] & ymax < B9[i]){ ylim2 <- c(0, B9[i]); break }
        if(ymax >= B9[i] & ymax < A1[i]){ ylim2 <- c(0, A1[i]*1.2); break }
      }  

      dev.new()
      layout( matrix(1:2, 1, 2, byrow=TRUE) )

      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )
      plot( x1, y1, main=main, cex.lab=1.5, cex.axis=1.5, 
            las=1, xaxs="i", yaxs="i", 
            xlab=xlab, ylab=ylab, 
            xlim=xlim1, ylim=ylim1, type="n" )
      polygon(x3, y3, col="grey70")
      points(x1, y1, cex=1.5, pch=16)
      lines(x3, y3, col=2) 
      lines(c(0, 1), c(0, 1), col=1)


      text(xloc, yloc, bquote(paste(hat(italic(c)), " = ", 
         .(sprintf("%.4f", round(par[1], 4))), sep="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-1.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(K)), {""}["1"], " = ", 
        .(sprintf("%.4f", round(par[2], 4))), sep="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-2.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(K)), {""}["2"], " = ", 
         .(sprintf("%.4f", round(par[3], 4))), sep="")), pos=4, cex=1.5, col=1)
      if(simpver==4){
        text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(a)), " = ", 
           .(sprintf("%.4f", round(par[4], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-4.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(b)), " = ", 
          .(sprintf("%.4f", round(par[5], 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-5.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
            .(sprintf("%.4f", round(r.sq, 4))), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-6.0*(ylim1[2]-ylim1[1])/10, 
           bquote( paste(italic(n), " = ", 
           .(length(x)), sep="") ), pos=4, cex=1.5, col=1)
      }
      if(simpver==5){
        text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
           .(round(r.sq, 4)), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-4.0*(ylim1[2]-ylim1[1])/10, 
           bquote( paste(italic(n), " = ", 
           .(length(x)), sep="") ), pos=4, cex=1.5, col=1)
      }

      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )

      plot(x, y, main=main, cex.lab=1.5, cex.axis=1.5,   
           las=1, xaxs="i", yaxs="i", xlab=expression(italic(x)), 
           ylab=expression(italic(y)), xlim=xlim2, ylim=ylim2, type="n")
      points(x, y, cex=1.5, pch=16)
      lines(x2, y2, col=2, lty=1, lwd=1) 

    }

    para.tab <- data.frame( Parameter = c(Names, "r.sq", "RSS", 
                    "sample.size", "GC"), Estimate = c(par, r.sq, RSS, 
                    length(x), GC) )
    if(par.list == "T" | par.list == "TRUE" | par.list == "True"){
        print(para.tab)
        cat("\n")
    }
    return(list( x1=x1, y1=y1, x=x, y=y, par=par, r.sq=r.sq, 
      RSS=RSS, sample.size=length(x), GC=GC) )

  }
}
