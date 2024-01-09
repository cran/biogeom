MPerformanceE <- function(P, x, simpver = 1){

  p <- length(P)

  if(is.null(simpver)){
    if(p != 7) 
    stop("The length of 'P' shoud equal 7 for the modified version!")
    c    <- P[1]
    K1   <- P[2]
    K2   <- P[3]
    xmin <- P[4]
    xmax <- P[5]
    a    <- P[6]
    b    <- P[7] 
  }

  if(!is.null(simpver)){

    if( !(simpver %in% seq(1, 5, by=1)) )
    stop("'simpver' should be chosen in versions 1 to 5!")  
  
    if(simpver==1){
      if(p != 6) 
      stop("The length of 'P' shoud equal 6 for the simplified version 1!")
      c    <- P[1]
      K1   <- P[2]
      K2   <- P[3]
      xmin <- 0
      xmax <- P[4]
      a    <- P[5]
      b    <- P[6] 
    }

   if(simpver==2){
      if(p != 5) 
      stop("The length of 'P' shoud equal 5 for the simplified version 2!")
      c    <- P[1]
      K1   <- P[2]
      K2   <- P[3]
      xmin <- P[4]
      xmax <- P[5]
      a    <- 1
      b    <- 1 
    }

    if(simpver==3){
      if(p != 4) 
      stop("The length of 'P' shoud equal 4 for the simplified version 3!")
      c    <- P[1]
      K1   <- P[2]
      K2   <- P[3]
      xmin <- 0
      xmax <- P[4]
      a    <- 1
      b    <- 1 
    }

    if(simpver==4){
      if(p != 5) 
      stop("The length of 'P' shoud equal 5 for the simplified version 4!")
      c    <- P[1]
      K1   <- P[2]
      K2   <- P[3]
      xmin <- 0
      xmax <- sqrt(2)
      a    <- P[4]
      b    <- P[5] 
    }

    if(simpver==5){
      if(p != 3) 
      stop("The length of 'P' shoud equal 3 for the simplified version 5!")
      c    <- P[1]
      K1   <- P[2]
      K2   <- P[3]
      xmin <- 0
      xmax <- sqrt(2)
      a    <- 1
      b    <- 1 
    }

  }

  fun0 <- function(y){
      y[y < xmin] <- xmin
      y[y > xmax] <- xmax
      return(y)
  }

  z <- fun0(x)

  temp1 <- c*( 1 - exp(-K1*(z-xmin)) )^a*( 1 - exp(K2*(z-xmax)) )^b 
  return(temp1)
}




