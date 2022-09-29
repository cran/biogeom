SurfaceAreaPE <- function(P, simpver = NULL, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
          stop.on.error = TRUE, keep.xy = FALSE, aux = NULL){

  if(P[1] < 0 | P[2] < 0)
    stop("a and b should be positive real numbers!")

  p <- length(P)
 
  if(is.null(simpver)){
    if(p != 5) 
      stop("The number of parameters is incorrect!")
  }

  if(!is.null(simpver)){

    if( !(simpver %in% seq(1, 3, by=1)) )
      stop("'simpver' should be chosen in versions 1 to 3!")
  
    if(simpver==1){
      if(p != 4) 
        stop("The number of parameters is incorrect!")
    }

    if(simpver==2){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
    }

    if(simpver==3){
      if(p != 3) 
        stop("The number of parameters is incorrect!")
    }
  }

  inner.fun <- function(x){
    2 * pi * EPE(P=P, x, simpver=simpver) * sqrt( 
      1 + DEPE(P=P, x, simpver=simpver)^2 )
  }

  integrate( inner.fun, -P[1], P[1], subdivisions = subdivisions,
             rel.tol = rel.tol, abs.tol = abs.tol,
             stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux )$value

}



