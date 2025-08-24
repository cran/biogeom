SurfaceAreaAM <- function(model = "Hybrid", P, upper = Inf, subdivisions = 100L,
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
 
    inner.fun <- function(x){
      2 * pi * x * sqrt(1+(alpha*beta*sinh(beta*x)+2*gamma*x)^2)
    }

    temp <- integrate( inner.fun, lower=0, upper=upper, subdivisions = subdivisions,
               rel.tol = rel.tol, abs.tol = abs.tol,
               stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value
    return( as.numeric( temp ) )
  }

  if(model == "Superparabola"){
    p <- length(P)
 
    if(p != 2) 
      stop("There are two model parameters for the superparabolic equation!")

    beta1 <- P[1]
    beta2 <- P[2]
 
    inner.fun <- function(x){
      term_sq <- (beta1 * beta2)^2 * abs(x)^(2 * (beta2 - 1))
      2 * pi * x * sqrt(1 + term_sq)
    }

    temp <- integrate( inner.fun, lower=0, upper=upper, subdivisions = subdivisions,
               rel.tol = rel.tol, abs.tol = abs.tol,
               stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value
    return( as.numeric( temp ) )
  }

}



