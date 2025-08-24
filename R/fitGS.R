
fitGS <- function(A, ini.val = NULL, control = list(), par.list = FALSE, fig.opt = TRUE){

  SA.temp <- A
  SA      <- sort( SA.temp )
  n       <- length(SA)
  x       <- 1:n
  y       <- SA
  z       <- cumsum(SA)  

  q1.est  <- mean(y[2:n]/y[1:(n-1)])  

  A1.est  <- 0 
  for(j in 1:n){
    A1.est <- A1.est + 1/n*y[j]/q1.est^(j-1)
  }

  obj.fun <- function(z){  
    A1 <- z[1]
    q  <- z[2]
    sum((A1*q^(x-1)-y)^2)
  }

  if(!is.null(ini.val)){
    ini.val <- as.list(ini.val)
    p <- length(ini.val)
    s <- 1
    for (i in 1:p) {
      s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat     <- matrix(NA, nrow = s, ncol = (p + 1))
  }
  if(is.null(ini.val)){
    ini.val <- t(cbind(c(A1.est, q1.est))) 
    p       <- 2
    mat     <- matrix(NA, nrow = 1, ncol = p+1)
  }

  for(i in 1:nrow(ini.val)){
    res      <- optim(ini.val[i, ], obj.fun, control = control)
    mat[i, ] <- c(res$par, res$val)
  }

  Names <- rep(NA, len = p)
  for (k in 1:p) {
      Names[k] <- paste("P[", k, "]", sep = "")
  }
  colnames(mat) <- c(Names, "RSS")
  ind <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
  MAT <- mat
  par <- as.vector(mat[ind, 1:p])
  PAR <- par[1:p]


  A.est   <- PAR[1]
  q.est   <- PAR[2]
  y.theo  <- A.est*q.est^(1:n-1) 
  z.theo  <- A.est*(1-q.est^(1:n))/(1-q.est)
  RSS     <- sum((y-y.theo)^2)
  r.sq    <- 1 - sum((y-y.theo)^2)/sum((y - mean(y))^2)
  RMSE1   <- sqrt( sum((y-y.theo)^2)/n )
  MAPE1   <- mean(abs(y-y.theo)/y * 100)
  RMSE2   <- sqrt( sum((z-z.theo)^2)/n )
  MAPE2   <- mean(abs(z-z.theo)/z * 100)

  para.tab <- data.frame(Parameter = c(Names, "r.sq", "RSS",  
      "sample.size", "RMSE1", "MAPE1", "RMSE2", "MAPE2"), 
      Estimate = c(par, r.sq, RSS, length(x), RMSE1, MAPE1, RMSE2, MAPE2))

  if(par.list == "T" | par.list == "TRUE" | par.list == "True"){
    print(para.tab)
    cat("\n")
  }

  if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True") {
 
    dev.new()
    layout( matrix(1:2, 1, 2, byrow=TRUE) )

    par( mar=c(5, 5, 2, 2) )
    par( mgp=c(3.5, 1, 0) )
    plot(x, y, cex.lab=1.5, cex.axis=1.5, cex=1.5, xlab=expression(italic(x)), 
         ylab=expression(paste("Sorted ", italic(y), " values in ascending order", sep="")))
    lines(x, y.theo, col=2)

    par( mar=c(5, 5, 2, 2) )
    par( mgp=c(3.5, 1, 0) )
    plot(x, z, cex.lab=1.5, cex.axis=1.5, cex=1.5, xlab=expression(italic(x)), 
         ylab=expression(paste("Cumulative sorted ", italic(y), " values", sep="")))
    lines(x, z.theo, col=2)

  }

  return(list( x=x, y=y, y.theo=y.theo, z=z, z.theo=z.theo, par=PAR, r.sq=r.sq, RSS=RSS, 
               sample.size=length(x), RMSE1=RMSE1, MAPE1=MAPE1, RMSE2=RMSE2, MAPE2=MAPE2 ))
}


