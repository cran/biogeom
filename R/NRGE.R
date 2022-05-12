
NRGE <- function(P, x){
    A     <- P[1]
    B     <- P[2]
    C     <- P[3]
    D     <- P[4]
    fun0 <- function(z){
        z[z < -A/2] <- -A/2
        z[z > A/2]  <- A/2
        return(z)
    }
    z1   <- fun0(x)
    temp <- B/2*sqrt((A^2-4*z1^2)/(A^2+8*C*z1+4*C^2))*
        (1-(sqrt(5.5*A^2+11*B*C+4*C^2)*
        (sqrt(3)*B*A-2*D*sqrt(A^2+2*C*A+4*C^2)))/
        (sqrt(3)*B*A*(sqrt(5.5*A^2+11*A*C+4*C^2)-2*
         sqrt(A^2+2*C*A+4*C^2)))*
         (1-sqrt((A*(A^2+8*C*z1+4*C^2))/(2*(A-2*C)*
         z1^2+(A^2+8*A*C-4*C^2)*z1+2*A*C^2+A^2*C+A^3)))) 
    return(temp)
}
