
areaGE <- function(expr, P, m = 1, simpver = NULL,  
            nval = 1, subdivisions = 100L,
            rel.tol = .Machine$double.eps^0.25, 
            abs.tol = rel.tol,
            stop.on.error = TRUE, keep.xy = FALSE, 
            aux = NULL){
  
   area.fun <- function(x){
     0.5*(expr(P, x, m=m, simpver=simpver, 
     nval=nval))^2
   } 

   temp <- integrate(area.fun, 0, 2*pi, 
           subdivisions = subdivisions,
           rel.tol = rel.tol, 
           abs.tol = abs.tol,
           stop.on.error = stop.on.error, 
           keep.xy = keep.xy, aux = aux)
   temp$value
}