EPE <- function(P, x, simpver = NULL){

  if(P[1] < 0 | P[2] < 0)
    stop("a and b should be positive real numbers!")
  if(min(x)[1] < -P[1] | max(x)[1] > P[1])
    stop("The value of x should be between -a and a!")

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

  d0 <- b/a
  d1 <- c1*(b/a)
  d2 <- c2*(b/a)
  d3 <- c3*(b/a)
  1/a^3*sqrt(a^2-x^2)*(a^3*d0 + a^2*d1*x + a*d2*x^2 + d3*x^3)

}



