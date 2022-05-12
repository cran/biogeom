adjdata <- function(x, y, ub.np = 2000, times = 1.2, 
               len.pro = 1/20, index.sp = 1){

  if(length(x)!=length(y)) 
    stop("The length of 'x' should equal that of 'y'!") 
  if(times < 1) 
    stop("The value of 'times' should be greater than 1!")
  if(len.pro > 1 | len.pro <= 0)
    stop("The value of 'len.pro' should range 0 to 1!")
  if( !(index.sp %in% 1:length(x)) )
    stop("'index.sp' should be an integer from 1 to length(x)!")

  start.mat  <- cbind(x,y)
  start.mat  <- na.omit(start.mat)
  x          <- start.mat[,1]
  y          <- start.mat[,2]

  ##### index.sp: the specified index of the starting point ####
  if(index.sp > 1){
    inde <- c(index.sp:length(x), 1:(index.sp-1))
    x    <- x[inde]
    y    <- y[inde]    
  }
  ##############################################################

  x.max  <- 2*max( abs(x) )[1]
  y.max  <- 2*max( abs(y) )[1]
  x      <- x + x.max
  y      <- y + y.max

  poly1  <- owin(xrange=c(min(x)[1], 
            max(x)[1]), yrange=c(min(y)[1],max(y)[1]))
  data0  <- ppp(x, y, window= poly1) 
  x      <- data0$x
  y      <- data0$y

  ini.np     <- ub.np*times
  if(length(x) > ini.np){
      ind   <- floor(as.numeric(quantile(1:length(x), 
                   probs=seq(0, 1, len=ini.np))))
      x     <- x[ind]
      y     <- y[ind]
      data0 <- ppp(x, y, window=poly1)   
  }

  dis.mat    <- pairdist(data0)
  len.max    <- max(dis.mat)[1]
  diag(dis.mat) <- NA 
  adj.mat    <- matrix(NA, nrow=data0$n, ncol=2)
  data1      <- cbind(x=data0$x, y=data0$y)
  dis.pool   <- c()
  for(i in 1:data0$n){
    if(i == 1){
      temp <- dis.mat[i, ]
      inde <- 1
    }
    old.inde     <- inde
    temp         <- dis.mat[inde, ]  
    if(length(temp[!is.na(temp)] ) > 0){ 
      inde <- which(temp == min(temp[!is.na(temp)]))[1]
      adj.mat[i, ] <- c( data1[old.inde, 1], data1[old.inde, 2] )
      dis.mat[old.inde, ] <- NA
      dis.mat[ ,old.inde] <- NA

      temp.dis <- min(temp[!is.na(temp)])[1]
      if(temp.dis > len.pro*len.max){
        adj.mat[i, ]    <- adj.mat[i-1, ]
        dis.mat[inde, ] <- NA
        dis.mat[, inde] <- NA
        inde            <- old.inde
      } 
    }
    if(length(temp[!is.na(temp)] ) == 0){ 
      adj.mat[i, ] <- c( data1[old.inde, 1], data1[old.inde, 2] )
      dis.mat[old.inde, ] <- NA
      dis.mat[ ,old.inde] <- NA
    }

  }
  adj.mat   <- unique(adj.mat)
  x         <- adj.mat[,1]
  y         <- adj.mat[,2] 
  poly.temp <- NULL
  try(poly.temp <- as.polygonal( owin(poly=list(x=x, 
                       y=y )) ), silent=TRUE)
  if( is.null(poly.temp) ){   
    x       <- x[length(x):1]
    y       <- y[length(y):1]
    poly0   <- as.polygonal( owin(poly=list(x=x, y=y )) )
  }
  if( !is.null(poly.temp) ){
    poly0   <- as.polygonal( owin(poly=list(x=x, y=y )) )
  }

  if(length(x) > ub.np){
      ind   <- floor(as.numeric(quantile(1:length(x), 
                   probs=seq(0, 1, len=ub.np))))
      x     <- x[ind]
      y     <- y[ind]
  }

  list(x=x-x.max, y=y-y.max)
}


