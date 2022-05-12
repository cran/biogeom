fitsigmoid <- function(expr, x, y, ini.val, simpver = 1, 
                control = list(), par.list = FALSE, fig.opt = FALSE,
                xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,         
                main = NULL, subdivisions = 100L,
                rel.tol=.Machine$double.eps^0.25, 
                abs.tol = rel.tol, stop.on.error = TRUE, 
                keep.xy = FALSE, aux = NULL){

    if(length(x)!=length(y)) 
        stop("'x' should have the same data length as 'y'!")
    Tem <- cbind(x, y)
    Tem <- na.omit(Tem)
    x   <- Tem[,1]
    y   <- Tem[,2]  

    if(!is.null(simpver)){
       if( !(simpver %in% seq(1, 3, by=1)) )
       stop("'simpver' should be chosen in versions 1 to 3!")    
    }

    ini.val <- as.list(ini.val)
    p       <- length(ini.val)
    s       <- 1
    for (i in 1:p) {
        s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat     <- matrix(NA, nrow = s, ncol = (p + 1))
  
    if( !identical(expr, MBriereE) ){ 
      obj.fun <- function(P){  
        if(is.null(simpver)){
          Lower <- P[3]
          Upper <- P[4]
          if(Lower >= Upper | P[1] <= 0 | P[2] <= Lower | P[2] >= Upper | P[5] <= 0)
            temp <- Inf
          if(Upper > Lower & P[1] > 0 & P[2] > Lower & P[2] < Upper & P[5] > 0){
            y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
            temp <- sum((y.theo-y)^2) 
          }
        }
        if(!is.null(simpver)){
          if(simpver==1){
            Lower <- 0
            Upper <- P[3]
            if(Lower >= Upper | P[1] <= 0 | P[2] <= Lower | P[2] >= Upper| P[4] <= 0)
                temp <- Inf
            if(Upper > Lower & P[1] > 0 & P[2] > Lower & P[2] < Upper & P[4] > 0){
              y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
              temp <- sum((y.theo-y)^2) 
            }
          }
          if(simpver==2){
            Lower <- P[3]
            Upper <- P[4]
            if(Lower >= Upper | P[1] <= 0 | P[2] <= Lower | P[2] >= Upper)
                temp <- Inf
            if(Upper > Lower & P[1] > 0 & P[2] > Lower & P[2] < Upper){
              y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
              temp <- sum((y.theo-y)^2) 
            }
          }
          if(simpver==3){
            Lower <- 0
            Upper <- P[3]
            if(Lower >= Upper | P[1] <= 0 | P[2] <= Lower | P[2] >= Upper)
                temp <- Inf
            if(Upper > Lower & P[1] > 0 & P[2] > Lower & P[2] < Upper){
              y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
              temp <- sum((y.theo-y)^2) 
            }
          }
        }
        return(temp) 
      }
    }

    if( identical(expr, MBriereE) ){
      obj.fun <- function(P){  
        if(is.null(simpver)){
          Lower <- P[3]
          Upper <- P[4]
          if(Lower >= Upper | Upper <= 0 | P[1] <= 0 | P[2] <= 0 | P[5] <= 0)
              temp <- Inf
          if(Upper > Lower & Upper > 0 & P[1] > 0 & P[2] > 0 & P[5] > 0){
            y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
            temp <- sum((y.theo-y)^2) 
          }
        }
        if(!is.null(simpver)){
          if(simpver==1){     
            Lower <- 0
            Upper <- P[3]
            if(Lower >= Upper | P[1] <= 0 | P[2] <= 0 | P[4] <= 0)
                temp <- Inf
            if(Upper > Lower & P[1] > 0 & P[2] > 0 & P[4] > 0){
              y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux)             
              temp <- sum((y.theo-y)^2)
            }
          }
          if(simpver==2){
            Lower <- P[3]
            Upper <- P[4]
            if(Lower >= Upper | Upper <= 0 | P[1] <= 0 | P[2] <= 0)
                temp <- Inf
            if(Upper > Lower & Upper > 0 & P[1] > 0 & P[2] > 0){
              y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
              temp <- sum((y.theo-y)^2) 
            }
          }
          if(simpver==3){
            Lower <- 0
            Upper <- P[3]
            if(Lower >= Upper | P[1] <= 0 | P[2] <= 0)
                temp <- Inf
            if(Upper > Lower & P[1] > 0 & P[2] > 0){
              y.theo <- sigmoid(expr, P, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, 
                abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux) 
              temp <- sum((y.theo-y)^2) 
            }
          }
        }
        return(temp) 
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
    ind  <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
    MAT  <- mat
    par  <- as.vector(mat[ind, 1:p])
    PAR  <- par






    y.theo <- sigmoid(expr, P=PAR, x, simpver=simpver,                
                subdivisions = subdivisions,
                rel.tol = rel.tol, abs.tol = abs.tol,
                stop.on.error = stop.on.error, 
                keep.xy = keep.xy, aux = aux)
    RSS     <- sum( (y - y.theo)^2 )
    r.sq    <- 1-sum((y-y.theo)^2)/sum((y-mean(y))^2)


    if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True"){  
        xinterv <- (max(x)[1]-min(x)[1])/10 
        xval    <- seq(min(x)[1]-xinterv, max(x)[1]+xinterv, len=2000)
        yval    <- sigmoid(expr, P=PAR, xval, simpver=simpver,                
                     subdivisions = subdivisions,
                     rel.tol = rel.tol, abs.tol = abs.tol,
                     stop.on.error = stop.on.error, 
                     keep.xy = keep.xy, aux = aux)
        if(is.null(xlim)){
          xlim     <- c(min(x)[1]-xinterv, max(x)[1]+xinterv)
        }

        if(!is.null(xlim)){
          if(length(xlim)!=2) stop("'xlim' should have two elements!")         
          xval    <- seq(xlim[1], xlim[2], len=2000)
          yval    <- sigmoid(expr, P=PAR, xval, simpver=simpver,                
                     subdivisions = subdivisions,
                     rel.tol = rel.tol, abs.tol = abs.tol,
                     stop.on.error = stop.on.error, 
                     keep.xy = keep.xy, aux = aux)
          xlim <- xlim
        }
        if(!is.null(ylim)){
          if(length(ylim)!=2) stop("'ylim' should have two elements!")
          ylim <- ylim
        }

        x_ran  <- c(x, xval)
        y_ran  <- c(y, yval) 

        if(is.null(xlab))
          xlab <- expression(italic(x)) 
        if(is.null(ylab))
          ylab <- expression(italic(y))      
      
        dev.new()
        plot( x_ran, y_ran, xlab=xlab, ylab=ylab, type="n", 
              xlim=xlim, ylim=ylim, cex.lab=1.5, cex.axis=1.5 )
        lines(xval, yval, type="l", asp=1, col=2, lwd=2)
        points( x, y, cex=1.5, col="grey40" )
        title(main=main, cex.main=1.5, col.main=4, font.main=1)
    } 

    para.tab <- data.frame( Parameter = c(Names, 
                    "r.sq", "RSS", "sample.size"), 
                    Estimate = c(par, r.sq, RSS, length(x)) )

    if(par.list == "T" | par.list == "TRUE" | par.list == "True" ){
        print(para.tab)
        cat("\n")
    }
    return(list( par=PAR, r.sq=r.sq, RSS=RSS, sample.size=length(x),
                 x=x, y=y, y.pred=y.theo ) )
}










