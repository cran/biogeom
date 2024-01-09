MbetaE <- function(P, x, simpver = 1){

  p <- length(P)

  if(is.null(simpver)){
    if(p != 5) 
    stop("The length of 'P' shoud equal 5 for the modified version!")
    yopt  <- P[1]
    xopt  <- P[2]
    xmin  <- P[3]
    xmax  <- P[4]
    delta <- P[5] 
  }

  if(!is.null(simpver)){

    if( !(simpver %in% seq(1, 3, by=1)) )
    stop("'simpver' should be chosen in versions 1 to 3!")  
  
    if(simpver==1){
      if(p != 4) 
      stop("The length of 'P' shoud equal 4 for the simplified version 1!")
      yopt  <- P[1]
      xopt  <- P[2]
      xmax  <- P[3]
      delta <- P[4]
      xmin  <- 0 
    }

   if(simpver==2){
      if(p != 4) 
      stop("The length of 'P' shoud equal 4 for the simplified version 2!")
      yopt  <- P[1]
      xopt  <- P[2]
      xmin  <- P[3]
      xmax  <- P[4]
      delta <- 1
    }

    if(simpver==3){
      if(p != 3) 
      stop("The length of 'P' shoud equal 3 for the simplified version 3!")
      yopt  <- P[1]
      xopt  <- P[2]
      xmax  <- P[3]
      xmin  <- 0
      delta <- 1
    }

  }

  fun0 <- function(y){
      y[y < xmin] <- xmin
      y[y > xmax] <- xmax
      return(y)
  }
  z <- fun0(x)

  yopt*(((xmax-z)/(xmax-xopt)*((z-xmin)/(xopt-xmin))^((xopt-xmin)/(xmax-xopt))))^delta

}









