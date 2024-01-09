
fitLorenz <- function(expr, z, ini.val, simpver = 4, 
               control = list(), par.list = FALSE, 
               fig.opt = FALSE, np = 2000, 
               xlab=NULL, ylab=NULL, main = NULL, subdivisions = 100L,
               rel.tol = .Machine$double.eps^0.25, 
               abs.tol = rel.tol, stop.on.error = TRUE, 
               keep.xy = FALSE, aux = NULL ){

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
      MPerformanceE(P=par, x, simpver=simpver)
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
      ylim2  <- c(0, 0.10)


      dev.new()
      layout( matrix(1:2, 1, 2, byrow=TRUE) )

      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )
      plot( x1, y1, main=main, cex.lab=1.5, cex.axis=1.5, 
            las=1, xaxs="i", yaxs="i", 
            xlab=xlab, ylab=ylab, 
            xlim=xlim1, ylim=ylim1, type="n" )
      polygon(x3, y3, col="grey70")
      lines(x3, y3, col=2) 
      lines(c(0, 1), c(0, 1), col=1)
      points(x1, y1, cex=1.5, pch=16)

      text(xloc, yloc, bquote(paste(hat(italic(c)), " = ", 
         .(round(par[1], 4)), sep="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-1.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(K)), {""}["1"], " = ", 
         .(round(par[2], 4)), sep="")), pos=4, cex=1.5, col=1)
      text(xloc, yloc-2.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(K)), {""}["2"], " = ", 
         .(round(par[3], 4)), sep="")), pos=4, cex=1.5, col=1)
      if(simpver==4){
        text(xloc, yloc-3.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(a)), " = ", 
           .(round(par[4], 4)), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-4.0*(ylim1[2]-ylim1[1])/10, bquote(paste(hat(italic(b)), " = ", 
           .(round(par[5], 4)), sep="")), pos=4, cex=1.5, col=1)
        text(xloc, yloc-5.0*(ylim1[2]-ylim1[1])/10, bquote(paste(italic(r), {""}^{"2"}, " = ", 
           .(round(r.sq, 4)), sep="")), pos=4, cex=1.5, col=1)
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
      lines(x2, y2, col=2, lty=1, lwd=1) 
      points(x, y, cex=1.5, pch=16)

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







