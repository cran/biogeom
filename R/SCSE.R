SCSE <- function(P, x){
    gamma <- P[1]
    alpha <- P[2]
    beta  <- P[3]
    x^gamma * (1-(1-x)^alpha)^beta
}
