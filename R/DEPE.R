DEPE <- function(P, x, simpver = NULL){

  if(P[1] < 0 | P[2] < 0)
    stop("a and b should be positive real numbers!")
  x[x < -P[1]] <- -P[1]
  x[x >  P[1]] <- P[1]
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
        c2       <- P[3]
        c3       <- 0
    }
  }

  -b/a^2*x*(1-(x/a)^2)^(-0.5)*(1+c1/a*x+c2*(x/a)^2+c3*(x/a)^3) +
    b*sqrt(1-(x/a)^2)*(c1/a+2*c2/a^2*x+3*c3/a^3*x^2)
}



