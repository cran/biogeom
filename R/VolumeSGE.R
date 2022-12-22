VolumeSGE <- function(P, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL){
  p <- length(P)
  if(p != 3) 
    stop("The number of parameters should equal 3!")
 
  inner.fun <- function(phi){    
    SGE <- GE(P=P, phi=phi, simpver=1, m=1)
    2/3*pi*SGE^3*sin(phi)
  }

  temp <- integrate( inner.fun, 0, pi, subdivisions = subdivisions,
             rel.tol = rel.tol, abs.tol = abs.tol,
             stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value
  as.numeric( temp )
}



