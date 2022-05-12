TGE <- function(P, phi, m = 1, simpver = NULL, nval = 1){
# Twin Gielis equation
  p <- length(P)
 
  if(is.null(simpver)){
    if(p != 6) 
    stop("The length of 'P' shoud equal 6 for the original TGE version!")
    a        <- P[1]
    b        <- P[2]
    c        <- P[3]
    k        <- P[4] 
    n2       <- P[5]
    n3       <- P[6] 
   }

  if(!is.null(simpver)){

    if( !(simpver %in% seq(1, 5, by=1)) )
  stop("'simpver' should be chosen in versions 1 to 5!")
  
    if(simpver==1){
      if(p != 4) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      b        <- P[2]
      c        <- P[3]
      k        <- 1
      n2       <- P[4]
      n3       <- n2  
    }

    if(simpver==2){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
      if(!is.numeric(nval))
        stop("'nval' should be a number!") 

      a        <- P[1]
      b        <- P[2]
      c        <- P[3]
      k        <- 1
      n2       <- nval
      n3       <- nval  
    }

    if(simpver==3){
      if(p != 5) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      b        <- P[2]
      c        <- P[3]      
      k        <- 1
      n2       <- P[4]
      n3       <- P[5]   
    }


    if(simpver==4){
      if(p != 5) 
        stop("The number of parameters is incorrect!")
      a        <- P[1]
      b        <- P[2]
      c        <- P[3]   
      k        <- P[4] 
      n2       <- P[5]
      n3       <- n2  
    }

    if(simpver==5){
      if(p != 4) 
        stop("The number of parameters is incorrect!")
      if(!is.numeric(nval))
        stop("'nval' should be a number!") 
      a        <- P[1]
      b        <- P[2]
      c        <- P[3]   
      k        <- P[4] 
      n2       <- nval
      n3       <- nval 
    }

  }

  re  <- (abs(cos(m/4*phi))^n2 + abs(1/k*sin(m/4*phi))^n3)^(-1) 
  r   <- exp( 1/(a+b*log(re)) + c )
  return( r )
}
