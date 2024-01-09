sigmoid <- function(expr, P, x, simpver = 1, subdivisions = 100L,
                    rel.tol = .Machine$double.eps^0.25, 
                    abs.tol = rel.tol, stop.on.error = TRUE, 
                    keep.xy = FALSE, aux = NULL){

   rate.fun <- function(x){
     expr(P, x=x, simpver=simpver)
   }

   if(is.null(simpver) & !identical(expr, MPerformanceE)){
     Lower <- P[3]
     Upper <- P[4]
   }

   if(is.null(simpver) & identical(expr, MPerformanceE)){
     Lower <- P[4]
     Upper <- P[5]
   }

   if(!is.null(simpver) & !identical(expr, MPerformanceE)){

    if( !(simpver %in% seq(1, 3, by=1)) )
    stop("'simpver' should be chosen in versions 1 to 3!")  
  
    if(simpver==1 | simpver==3){
      Lower <- 0
      Upper <- P[3]
    }
    if(simpver==2){
      Lower <- P[3]
      Upper <- P[4]
    }
   }

   if(!is.null(simpver) & identical(expr, MPerformanceE)){

    if( !(simpver %in% seq(1, 5, by=1)) )
    stop("'simpver' should be chosen in versions 1 to 5!")  
  
    if(simpver==1 | simpver==3){
      Lower <- 0
      Upper <- P[4]
    }
    if(simpver==2){
      Lower <- P[4]
      Upper <- P[5]
    }
    if(simpver==4 | simpver==5){
      Lower <- 0
      Upper <- sqrt(2)
    }
   }

   yval <- c()
   for(i in 1:length(x)){
     if(x[i] < Lower) x[i] <- Lower
     if(x[i] > Upper) x[i] <- Upper
     temp <- integrate(rate.fun, lower=Lower, upper=x[i], 
               subdivisions = subdivisions,
               rel.tol = rel.tol, 
               abs.tol = abs.tol,
               stop.on.error = stop.on.error, 
               keep.xy = keep.xy, aux = aux)
     yval <- c(yval, temp$value)  
   }
  
   return( yval )

}
