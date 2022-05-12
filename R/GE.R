GE <- function(P, phi, m = 1, simpver = NULL, nval = 1){

  p <- length(P)
 
  if(is.null(simpver)){
    if(p != 5) 
      stop("The number of parameters is incorrect!")
    a        <- P[1]
    k        <- P[2]
    n1       <- P[3]
    n2       <- P[4]
    n3       <- P[5]   
  }

  if(!is.null(simpver)){

    if( !(simpver %in% seq(1, 9, by=1)) )
  stop("'simpver' should be chosen in versions 1 to 9!")
  
    if(simpver==1){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      k        <- 1
      n1       <- P[2]
      n2       <- P[3]
      n3       <- n2  
    }

    if(simpver==2){
      if(p != 2) 
        stop("The number of parameters is incorrect!")
    if(!is.numeric(nval))
      stop("'nval' should be a number!") 
      a        <- P[1]
      k        <- 1
      n1       <- P[2]
      n2       <- nval
      n3       <- nval  
    }


    if(simpver==3){
      if(p != 1) 
        stop("The number of parameters is incorrect!")
      if(!is.numeric(nval))
        stop("'nval' should be a number!") 
      a        <- P[1]
      k        <- 1
      n1       <- nval
      n2       <- nval
      n3       <- nval  
    }

    if(simpver==4){
      if(p != 2) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      k        <- 1
      n1       <- P[2]
      n2       <- n1
      n3       <- n2   
    }

    if(simpver==5){
      if(p != 4) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      k        <- 1 
      n1       <- P[2]
      n2       <- P[3]
      n3       <- P[4]  
    }

    if(simpver==6){
      if(p != 4) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      k        <- P[2] 
      n1       <- P[3]
      n2       <- P[4]
      n3       <- n2  
    }

    if(simpver==7){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
      if(!is.numeric(nval))
        stop("'nval' should be a number!") 
      a        <- P[1]
      k        <- P[2] 
      n1       <- P[3]
      n2       <- nval
      n3       <- nval 
    }

    if(simpver==8){
      if(p != 2) 
        stop("The number of parameters is incorrect!")
      if(!is.numeric(nval))
        stop("'nval' should be a number!") 
      a        <- P[1]
      k        <- P[2] 
      n1       <- nval
      n2       <- nval
      n3       <- nval
    }


    if(simpver==9){
      if(p != 3) 
  stop("The length of 'P' shoud equal 3 for the simplified version 9!")
      a        <- P[1]
      k        <- P[2] 
      n1       <- P[3]
      n2       <- n1
      n3       <- n2
    }

  }

  r  <- a*(abs(cos(m/4*phi))^n2 + abs(1/k*sin(m/4*phi))^n3)^(-1/n1) 
  return( r )
}
