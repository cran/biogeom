TSE <- function(P, x, simpver=NULL){
    temp       <- 1-x^2
    ind0       <- which(temp < 0)
    temp[ind0] <- 0
    z0         <- sqrt(temp)
    z1         <- x*z0
    z2         <- x^2*z0
    z3         <- x^3*z0
    if(is.null(simpver)){
      d0 <- P[1]
      d1 <- P[2]
      d2 <- P[3]
      d3 <- P[4]
      y  <- d0*z0 + d1*z1 + d2*z2 + d3*z3   
    }
    if(!is.null(simpver) && !is.numeric(simpver))
      stop("'simpver' must be chosen from 1, 2 and 3")  
    if(!is.null(simpver) && simpver==1){    
      d0 <- P[1]
      d1 <- P[2]
      d2 <- P[3]
      y <- d0*z0 + d1*z1 + d2*z2       
    }
    if(!is.null(simpver) && simpver==2){    
      d0 <- P[1]
      d1 <- P[2]
      y <- d0*z0 + d1*z1 
    }
    if(!is.null(simpver) && simpver==3){   
      d0 <- P[1]
      d2 <- P[2] 
      y <- d0*z0 + d2*z2
    }
    return(y)
}




