
fitovate <- function(expr, x, y, ini.val, 
              par.list = FALSE, stand.fig = TRUE, control = list(), 
              angle = NULL, fig.opt = FALSE, index.xmax = 3, np = 2000, 
              xlim = NULL, ylim = NULL, unit = NULL, main = NULL){

    if(length(x)!=length(y)) 
        stop("'x' should have the same data length as 'y'!")
    Tem <- cbind(x, y)
    Tem <- na.omit(Tem)
    x   <- Tem[,1]
    y   <- Tem[,2]    

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
        yy1     <- y.obs[cond1] 
        xx2     <- x.obs[cond2]
        yy2     <- y.obs[cond2]  
 
        ind1    <- sort(xx1, decreasing=TRUE, index.return=TRUE)$ix
        xx3     <- xx1[ind1]
        yy3     <- yy1[ind1] 
        ind2    <- sort(xx2, decreasing=FALSE, index.return=TRUE)$ix  
        xx4     <- xx2[ind2]
        yy4     <- yy2[ind2]  

        res1    <- UpperFun(z[4:p], x=xx3)
        res2    <- LowerFun(z[4:p], x=xx4) 
        x.theo  <- c(res1$x, res2$x)
        y.theo  <- c(res1$y, res2$y)
        x.obs   <- c(xx3, xx4)
        y.obs   <- c(yy3, yy4)
        Value   <- sum( (y.obs - y.theo)^2 )
        #### Diminishes the abnormal values ###################
        if(identical(expr, MbetaE) | identical(expr, 
            MLRFE) | identical(expr, MBriereE) ){
            A <- z[4]
            B <- z[5]
            C <- z[6]
            if(A <= 0 | B < 0 | C <= 0) 
            Value <- Inf 
          }

        if(TRUE){
          if( identical(expr, MPerformanceE) ){
            A <- z[4] # c
            B <- z[5] # K1
            C <- z[6] # K2
            D <- z[7] # xmax
            if(A <= 0 | B < 0 | C < 0 | D <= 0) 
            Value <- Inf 
          }
        }   
        ######################################################        
        return( Value )       
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
    MAT  <- mat
    par  <- as.vector(mat[ind, 1:p])
    PAR  <- par[4:p]

    #### To solve the issue of abs(theta) > 2*pi ####
    if(par[3] > 2*pi){
        par[3] <- par[3]%%(2*pi)
    }
    if(par[3] < -2*pi){
        par[3] <- -((-par[3])%%(2*pi))
    }
    ##################################################

    goal.x0     <- par[1]
    goal.y0     <- par[2]
    goal.theta  <- par[3]  
    x.obs       <- x - goal.x0
    y.obs       <- y - goal.y0
    x.new       <- x.obs*cos(goal.theta) + y.obs*sin(goal.theta)
    y.new       <- y.obs*cos(goal.theta) - x.obs*sin(goal.theta)
    x.obs       <- x.new
    y.obs       <- y.new
    cond1       <- y.obs >= 0
    cond2       <- y.obs <  0
    xx1         <- x.obs[cond1]
    xx2         <- x.obs[cond2]
    yy1         <- y.obs[cond1] 
    yy2         <- y.obs[cond2]  
    ind1        <- sort(xx1, decreasing=TRUE, index.return=TRUE)$ix
    xx3         <- xx1[ind1]
    yy3         <- yy1[ind1] 
    ind2        <- sort(xx2, decreasing=FALSE, index.return=TRUE)$ix  
    xx4         <- xx2[ind2]
    yy4         <- yy2[ind2]  
    x1          <- c(xx3, xx4)
    y1          <- c(yy3, yy4)
    res1        <- UpperFun(PAR, x=xx3)
    res2        <- LowerFun(PAR, x=xx4)
    x.theo      <- c(res1$x, res2$x)
    y.theo      <- c(res1$y, res2$y)
    x.obs       <- c(xx3, xx4)
    y.obs       <- c(yy3, yy4)
    x2          <- x.theo
    y2          <- y.theo
    r.obs       <- sqrt( x1^2 + y1^2 )
    phi.obs     <- acos(x1/r.obs)
    phiA        <- phi.obs 
    x.arr       <- seq(0, PAR[index.xmax], len=np)
    if( is.null(angle) ){        
        x.new    <- x1*cos(goal.theta) - y1*sin(goal.theta) + par[1]
        y.new    <- y1*cos(goal.theta) + x1*sin(goal.theta) + par[2] 
        x.theory <- x2*cos(goal.theta) - y2*sin(goal.theta) + par[1]
        y.theory <- y2*cos(goal.theta) + x2*sin(goal.theta) + par[2] 
        results2 <- curveovate(expr, P=par, x=x.arr)
        phiB     <- phiA + goal.theta        
    }
    if( !is.null(angle) ){        
        par.new    <- par  
        #### The third element of 'par.new' is compulsorily defined as 'angle' ####
        par.new[3] <- angle 
        ###########################################################################
        x.new    <- x1*cos(angle) - y1*sin(angle) + par.new[1]
        y.new    <- y1*cos(angle) + x1*sin(angle) + par.new[2]         
        x.theory <- x2*cos(angle) - y2*sin(angle) + par[1]
        y.theory <- y2*cos(angle) + x2*sin(angle) + par[2] 
        results2 <- curveovate(expr, P=par.new, x=x.arr)
        phiB     <- phiA + angle 
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
          plot(x1, y1, xlab=xlabel, ylab=ylabel, las=1,
               cex.lab=1.5, cex.axis=1.5, 
               type="l", lwd=3, col="grey50", asp=1)
          lines(x2, y2, col=2, lwd=2)
          title(main=main, cex.main=1.5, col.main=4, font.main=1)
          abline(h=0, lty=2, col=4)
          abline(v=0, lty=2, col=4)
    }

    if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True"){  
      x_ran  <- c(x.new, x.theory)
      y_ran  <- c(y.new, y.theory)      
      if(is.null(angle)){ 
          dev.new()
          plot( x_ran, y_ran, xlab=xlabel, ylab=ylabel, type="n", 
                asp=1, xlim=xlim, ylim=ylim, cex.lab=1.5, cex.axis=1.5 )
          points( x.new, y.new, cex=0.5, col="grey40" )
          lines(results2$x, results2$y, type="l", asp=1, col=2, lwd=2)
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
          plot( x_ran, y_ran, xlab=xlabel, ylab=ylabel, type="n", 
                asp=1, xlim=xlim, ylim=ylim, cex.lab=1.5, cex.axis=1.5 )
          points( x.new, y.new, cex=0.5, col="grey40" )
          lines(results2$x, results2$y, type="l", asp=1, col=2, lwd=2)
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

    L1      <- max(x1)[1] - min(x1)[1]
    W1      <- max(y1)[1] - min(y1)[1]  
    poly1   <- as.polygonal( owin(poly=list(x=x1, y=y1)) )
    P1      <- perimeter(poly1)
    A1      <- area(poly1)
    RES1    <- UpperFun(PAR, x=seq(0, PAR[3], len=np))
    RES2    <- LowerFun(PAR, x=seq(0, PAR[3], len=np))
    x3      <- c(RES1$x, RES2$x)
    y3      <- c(RES1$y, RES2$y)
    L2      <- max(x3)[1] - min(x3)[1]
    W2      <- max(y3)[1] - min(y3)[1] 
    poly2   <- as.polygonal( owin(poly=list(x=x3, y=y3)) )
    P2      <- perimeter(poly2)
    A2      <- area(poly2)
    RSS     <- sum((y.obs - y.theo)^2)
    r.sq    <- 1-sum((y.obs-y.theo)^2)/sum((y.obs-mean(y.obs))^2)

    para.tab <- data.frame( Parameter = c(Names, 
                    "r.sq", "RSS", "sample.size"), 
                    Estimate = c(par, r.sq, RSS, length(x)) )
    if(par.list == "T" | par.list == "TRUE"){
        print(para.tab)
        cat("\n")
    }
    return(list( par=par, r.sq=r.sq, RSS=RSS, sample.size=length(x),
                 scan.length=L1, scan.width=W1, scan.perimeter=P1, scan.area=A1,
                 pred.length=L2, pred.width=W2, pred.perimeter=P2, pred.area=A2,
                 x.stand.obs=x1, x.stand.pred=x2, 
                 y.stand.obs=y1, y.stand.pred=y2, 
                 x.obs=x.new, x.pred=x.theory, y.obs=y.new, y.pred=y.theory))
}
