PE <- function(P, zeta, simpver = NULL){

  p <- length(P)
 
  if(is.null(simpver)){
    if(p != 5) 
      stop("The number of parameters is incorrect!")
    a        <- P[1]
    b        <- P[2]
    c1       <- P[3]
    c2       <- P[4]
    c3       <- P[5]
  }

  if(!is.null(simpver)){

    if( !(simpver %in% seq(1, 3, by=1)) )
      stop("'simpver' should be chosen in versions 1 to 3!")
  
    if(simpver==1){
      if(p != 4) 
        stop("The number of parameters is incorrect!")
        a        <- P[1]
        b        <- P[2]
        c1       <- P[3]
        c2       <- P[4]
        c3       <- 0
    }

    if(simpver==2){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
        a        <- P[1]
        b        <- P[2]
        c1       <- P[3]
        c2       <- 0
        c3       <- 0
    }

    if(simpver==3){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
        a        <- P[1]
        b        <- P[2]
        c1       <- 0
        c2       <- P[4]
        c3       <- 0
    }
  }

  x <- b*cos(zeta)*(1+c1*sin(zeta)+c2*(sin(zeta))^2+c3*(sin(zeta))^3)
  y <- a*sin(zeta)
  r <- sqrt(x^2 + y^2)
  list( x=x, y=y, r=r )
}



