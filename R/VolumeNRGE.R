VolumeNRGE <- function(P, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL){

  if(P[1] <= 0)
    stop("A should be positive real numbers!")

  p <- length(P)
  if(p != 4) 
    stop("The number of parameters should equal 4!")

  inner.fun <- function(x){
    pi * NRGE(P=P, x)^2 
  }

  temp <- integrate( inner.fun, -P[1]/2, P[1]/2, subdivisions = subdivisions,
             rel.tol = rel.tol, abs.tol = abs.tol,
             stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value
  as.numeric( temp )

}



