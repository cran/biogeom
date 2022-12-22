curveETE <- function (P, np = 5000, 
    fig.opt = FALSE, deform.fun = NULL, Par = NULL, 
    xlim = NULL, ylim = NULL, unit = NULL, main = "") 
{
    if ((is.null(deform.fun) & !is.null(Par)) | (!is.null(deform.fun) & 
        is.null(Par))) 
        stop("'Par' should be provided when 'deform.fun' is not null.")

    x0    <- P[1]
    y0    <- P[2]
    theta <- P[3]
    npar  <- length(P)
    if(np %% 2 == 1) np <- np + 1
    xU    <- seq(P[4], -P[4], len = np/2)
    yU    <- ETE(P = P[4:npar], x = xU)
    xL    <- seq(-P[4], P[4], len = np/2)
    yL    <- -ETE(P = P[4:npar], x = xL)
    x     <- c(xU, xL)
    y     <- c(yU, yL)

    if (!is.null(deform.fun)) {
        Resu <- deform.fun(Par = Par, z = cbind(x, y))
        x <- Resu$x
        y <- Resu$y
    }
    x.rot    <- x * cos(theta) - y * sin(theta)
    y.rot    <- y * cos(theta) + x * sin(theta)
    r.rot    <- sqrt(x.rot^2 + y.rot^2)
    x.coordi <- x.rot + x0
    y.coordi <- y.rot + y0
    if (fig.opt == "T" | fig.opt == "TRUE") {
        if (is.null(xlim)) 
            xlim <- NULL
        if (is.null(ylim)) 
            ylim <- NULL
        if (!is.null(xlim)) 
            xlim <- xlim
        if (!is.null(ylim)) 
            ylim <- ylim
        if (!is.null(unit)) {
            xlabel <- bquote(paste(italic("x"), " (", .(unit), 
                ")", sep = ""))
            ylabel <- bquote(paste(italic("y"), " (", .(unit), 
                ")", sep = ""))
        }
        if (is.null(unit)) {
            xlabel <- bquote(italic("x"))
            ylabel <- bquote(italic("y"))
        }
        dev.new()
        plot(x.coordi, y.coordi, asp = 1, xlab = xlabel, ylab = ylabel, 
            pch = 1, cex = 2, cex.lab = 1.5, cex.axis = 1.5, 
            type = "l", lwd = 2, xlim = xlim, ylim = ylim)
        abline(h = y0, lty = 2, col = 4)
        abline(v = x0, lty = 2, col = 4)
        title(main = main, cex.main = 1.5, col.main = 4, font.main = 1)
        xu <- P[4]
        yu <- ETE(P = P[4:npar], x = xu)
        if (!is.null(deform.fun)) {
            Right <- deform.fun(Par = Par, z = cbind(xu, yu))
            xu <- Right$x
            yu <- Right$y
        }
        xv <- xu * cos(theta) - yu * sin(theta)
        yv <- yu * cos(theta) + xu * sin(theta)
        xv <- xv + x0
        yv <- yv + y0
        temp <- xv - x0
        if (temp != 0) {
            slope <- (yv - y0)/(xv - x0)
            abline(-slope * x0 + y0, slope, col = 1, lty = 2)
        }
        if (temp == 0) {
            abline(v = x0, col = 1, lty = 2)
        }
        points(x0, y0, cex = 2, pch = 16, col = 2)
        points(xv, yv, cex = 2, pch = 1, col = 4)
    }
    list(x = x.coordi, y = y.coordi)
}
