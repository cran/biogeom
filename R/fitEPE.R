
fitEPE <- function(x, y, ini.val, simpver = NULL, 
           control = list(), par.list = FALSE, 
           stand.fig = TRUE, angle = NULL, fig.opt = FALSE, np = 2000,
           xlim = NULL, ylim = NULL, unit = NULL, main = NULL ){
    if(length(x)!=length(y)) 
        stop("'x' should have the same data length as 'y'!")
    if(np %% 2 == 1) np <- np + 1
    Tem   <- cbind(x, y)
    Tem   <- na.omit(Tem)
    x     <- Tem[,1]
    y     <- Tem[,2] 
    ini.val <- as.list(ini.val)
    p       <- length(ini.val)
    s       <- 1
    for (i in 1:p) {
        s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat     <- matrix(NA, nrow = s, ncol = (p + 1))
    obj.fun <- function(z){
        x0      <- z[1]
        y0      <- z[2]
        theta   <- z[3]
        x.obs   <- x - x0
        y.obs   <- y - y0
        x.temp  <- x.obs*cos(theta) + y.obs*sin(theta)
        y.temp  <- y.obs*cos(theta) - x.obs*sin(theta)
        x.obs   <- x.temp
        y.obs   <- y.temp
        cond1   <- y.obs >= 0
        cond2   <- y.obs <  0
        xx1     <- x.obs[cond1]
        xx2     <- x.obs[cond2]     
        yy1     <- y.obs[cond1] 
        yy2     <- y.obs[cond2] 
  
        if(z[4] <= 0 | z[5] <= 0) RSS <- Inf 
        if(z[4] > 0 & z[5] > 0){
          if(length(xx1)==0 | length(xx2)==0) 
            RSS <- Inf
          if(length(xx1) > 0 & length(xx2) > 0){ 
            yU  <- EPE(P=z[4:p], x=xx1, simpver=simpver) 
            yL  <- -EPE(P=z[4:p], x=xx2, simpver=simpver) 
            RSS <- sum( (yU - yy1)^2 ) + sum( (yL - yy2)^2 )
          }
        }
        return( RSS )       
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

    ind  <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
    par  <- as.vector(mat[ind, 1:p])
    #### To solve the issue of abs(theta) > 2*pi ####
    if(par[3] > 2*pi){
        par[3] <- par[3]%%(2*pi)
    }
    if(par[3] < -2*pi){
        par[3] <- -((-par[3])%%(2*pi))
    }
    ##################################################
    goal.x0    <- par[1]
    goal.y0    <- par[2]
    goal.theta <- par[3]
    x.obs <- x - goal.x0
    y.obs <- y - goal.y0
    x.new <- x.obs*cos(goal.theta) + y.obs*sin(goal.theta)
    y.new <- y.obs*cos(goal.theta) - x.obs*sin(goal.theta)
    x.obs <- x.new
    y.obs <- y.new
    r.obs   <- sqrt( x.obs^2 + y.obs^2 )
    phi.obs <- acos(x.obs/r.obs) 
    cond1 <- y.obs >= 0
    cond2 <- y.obs <  0
    xx1   <- x.obs[cond1]
    xx2   <- x.obs[cond2]
    yy1   <- y.obs[cond1] 
    yy2   <- y.obs[cond2]
    yU      <- EPE(P=par[4:p], x=xx1, simpver=simpver) 
    yL      <- -EPE(P=par[4:p], x=xx2, simpver=simpver)       
    phi1    <- phi.obs[cond1]  
    phi2    <- 2*pi - phi.obs[cond2]
    phi.obs <- c(phi1, phi2)
    x.obs   <- c(xx1, xx2)
    y.obs   <- c(yy1, yy2)
    y2      <- c(yU, yL) 
    resu    <- sort(phi.obs, decreasing = FALSE, index.return= T )
    phi.obs <- resu$x
    Index   <- resu$ix
    x1      <- x.obs[Index]
    y1      <- y.obs[Index]
    y2      <- y2[Index]
    CombTab <- cbind(x1, y1, y2)
    CombTab <- na.omit(CombTab)
    x1      <- CombTab[,1]
    y1      <- CombTab[,2]
    y2      <- CombTab[,3]     

    xv1    <- seq(par[4], -par[4], len=np/2)
    xv2    <- seq(-par[4], par[4], len=np/2)
    yU.val <- EPE(P=par[4:p], x=xv1, simpver=simpver)
    yL.val <- -EPE(P=par[4:p], x=xv2, simpver=simpver)
    xv3 <- c(xv1, xv2)
    yv3 <- c(yU.val, yL.val)
    if( is.null(angle) ){ 
        x.new  <- x1*cos(goal.theta) - y1*sin(goal.theta) + par[1]
        y.new  <- y1*cos(goal.theta) + x1*sin(goal.theta) + par[2] 
        xv4    <- xv3*cos(goal.theta) - yv3*sin(goal.theta) + par[1]
        yv4    <- yv3*cos(goal.theta) + xv3*sin(goal.theta) + par[2] 
        yv5    <- y2*cos(goal.theta) + x1*sin(goal.theta) + par[2] 
    }
    if( !is.null(angle) ){        
        par.new    <- par  
        #### The third element of 'par.new' is compulsorily defined as 'angle' ####
        par.new[3] <- angle 
        ###########################################################################
        x.new    <- x1*cos(angle) - y1*sin(angle) + par.new[1]
        y.new    <- y1*cos(angle) + x1*sin(angle) + par.new[2] 
        xv4      <- xv3*cos(angle) - yv3*sin(angle) + par.new[1]
        yv4      <- yv3*cos(angle) + xv3*sin(angle) + par.new[2] 
        yv5      <- y2*cos(angle) + x1*sin(angle) + par.new[2] 
    }
    if(is.null(xlim)) xlim <- NULL
    if(is.null(ylim)) ylim <- NULL
    if(!is.null(xlim)) xlim <- xlim
    if(!is.null(ylim)) ylim <- ylim        
    if(!is.null(unit)){
      xlabel <- bquote( paste(italic("x"), " (", .(unit), ")", sep="") ) 
      ylabel <- bquote( paste(italic("y"), " (", .(unit), ")", sep="") )
    }
    if(is.null(unit)){
      xlabel <- bquote( italic("x") ) 
      ylabel <- bquote( italic("y") )
    }
    if(stand.fig == "T" | stand.fig == "TRUE" | stand.fig == "True"){
          dev.new()
          plot(x1, y1, xlab=xlabel, ylab=ylabel, 
               cex.lab=1.5, cex.axis=1.5, 
               type="l", lwd=3, col="grey50", asp=1)
          lines(x1, y2, col=2, lwd=2)
          title(main=main, cex.main=1.5, col.main=4, font.main=1)
          abline(h=0, lty=2, col=4)
          abline(v=0, lty=2, col=4)
    }
    if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True"){
      if(is.null(angle)){
            dev.new()
            plot( x.new, y.new, xlab=xlabel, ylab=ylabel, type="l",
                  cex.lab=1.5, cex.axis=1.5, col="grey50", lwd=3, asp=1, 
                  xlim=xlim, ylim=ylim )
            lines(xv4, yv4, type="l", asp=1, col=2, lwd=2)
            title(main=main, cex.main=1.5, col.main=4, font.main=1)
            if(abs(par[3])%%pi/2 != 0){
              slope <- tan(par[3])
              abline(-slope*par[1]+par[2], slope, col=1, lty=2)
              abline(v = par[1], col = 4, lty=2, lwd=1)
              abline(h = par[2], col = 4, lty=2, lwd=1)

            }
            if(abs(par[3])%%pi/2 == 0){
              abline(v = par[1], col = 1, lty=2)
              abline(v = par[1], col = 4, lty=2, lwd=1)
              abline(h = par[2], col = 4, lty=2, lwd=1)
            }
      } 
      if(!is.null(angle)){
            dev.new()
            plot( x.new, y.new, xlab=xlabel, ylab=ylabel, type="l", 
                  cex.lab=1.5, cex.axis=1.5, col="grey50", lwd=3, asp=1,
                  xlim=xlim, ylim=ylim )
            lines(xv4, yv4, type="l", asp=1, col=2, lwd=2)
            title(main=main, cex.main=1.5, col.main=4, font.main=1)
            if(abs(par.new[3])%%pi/2 != 0){
              slope <- tan(par.new[3])
              abline(-slope*par[1]+par[2], slope, col=1, lty=2)
              abline(v = par[1], col = 4, lty=2, lwd=1)
              abline(h = par[2], col = 4, lty=2, lwd=1)

            }
            if(abs(par.new[3])%%pi/2 == 0){
              abline(v = par[1], col = 1, lty=2)
              abline(v = par[1], col = 4, lty=2, lwd=1)
              abline(h = par[2], col = 4, lty=2, lwd=1)
            }
      }
    }
    x.range <- range(x1)
    y.range <- range(y1)
    length0 <- x.range[2] - x.range[1]
    width0  <- y.range[2] - y.range[1]    
    poly0   <- as.polygonal(owin(poly = list(x = x1, y = y1)))
    area0   <- area.owin(poly0)
    perimeter0 <- perimeter(poly0)
    RSS      <- sum((y1 - y2)^2)
    r.sq     <- 1-sum((y1-y2)^2)/sum((y1-mean(y1))^2)
    para.tab <- data.frame( Parameter = c(Names, 
                    "Length", "Width", "Area", "Perimeter", "r.sq", "RSS", 
                    "sample.size"), Estimate = c(par, length0, width0, 
                    area0, perimeter0, r.sq, RSS, length(x1)) )
    if(par.list == "T" | par.list == "TRUE" | par.list == "True"){
        print(para.tab)
        cat("\n")
    }
    return(list( par=par, scan.length=length0, scan.width=width0, 
      scan.area=area0, scan.perimeter=perimeter0, r.sq=r.sq, 
      RSS=RSS, sample.size=length(x1),      
      x.stand.obs=x1, y.stand.obs=y1, y.stand.pred=y2, 
      x.obs=x.new, y.obs=y.new, y.pred=yv5) )
}







