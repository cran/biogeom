VolumeAM <- function(model = "Hybrid", P, upper = Inf, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL){


  if(!(model %in% c("Hybrid", "Superparabola"))){
    stop("'model' only can be 'Hybrid' or 'Superparabola' here!")
  }

  if(model == "Hybrid"){
    p <- length(P)
 
    if(p != 3) 
      stop("There are three model parameters for the hybrid catenary-parabolic equation!")

    alpha <- P[1]
    beta  <- P[2]
    gamma <- P[3]
 
    Hybrid <- function(x){
      alpha * cosh(beta*x) + gamma*x^2 - alpha
    } 

    inner.fun <- function(x){
      2 * pi * abs(x * Hybrid(x))     
    }

    temp <- integrate( inner.fun, lower=0, upper=upper, subdivisions = subdivisions,
               rel.tol = rel.tol, abs.tol = abs.tol,
               stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value

    if( Hybrid(upper) > Hybrid(0) )
        return( pi * upper^2 * Hybrid(upper) - as.numeric( temp ) )
    if( Hybrid(upper) < Hybrid(0) )
        return( as.numeric( temp ) )
    if( Hybrid(upper) == Hybrid(0) )
        return(0)
  }

  if(model == "Superparabola"){
    p <- length(P)
 
    if(p != 2) 
      stop("There are two model parameters for the superparabolic equation!")

    beta1 <- P[1]
    beta2 <- P[2]

    Superparabola <- function(x){
      beta1 * abs(x)^beta2
    }
 
    inner.fun <- function(x){
      2 * pi * abs(x * Superparabola(x))      
    }

    temp <- integrate( inner.fun, lower=0, upper=upper, subdivisions = subdivisions,
               rel.tol = rel.tol, abs.tol = abs.tol,
               stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value

    

    if( Superparabola(upper) > Superparabola(0) )
        return( pi * upper^2 * Superparabola(upper) - as.numeric( temp ) )
    if( Superparabola(upper) < Superparabola(0) )
        return( as.numeric( temp ) )
    if( Superparabola(upper) == Superparabola(0) )
        return(0)
  }

}



