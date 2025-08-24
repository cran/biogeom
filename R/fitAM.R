
fitAM <- function (expr, x, y, ini.val, method = "Nelder-Mead", control = list(), 
    lower = -Inf, upper = Inf, par.list = FALSE, stand.fig = TRUE, fig.opt = FALSE,  
    angle = NULL, xlim = NULL, ylim = NULL, unit = NULL, main = NULL){

    if( !( method %in% c("Nelder-Mead", "L-BFGS-B") ) ){
        stop("'method' is only selected as Nelder-Mead or L-BFGS-B!")
    } 

    if (length(x) != length(y)) 
        stop("'x' should have the same data length as 'y'!")
    Tem <- cbind(x, y)
    Tem <- na.omit(Tem)
    x   <- Tem[, 1]
    y   <- Tem[, 2]
    ini.val <- as.list(ini.val)
    p <- length(ini.val)
    s <- 1
    for (i in 1:p) {
        s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat <- matrix(NA, nrow = s, ncol = (p + 1))
    obj.fun <- function(z) {
        x0     <- z[1]
        y0     <- z[2]
        theta  <- z[3]
        x.obs  <- x - x0
        y.obs  <- y - y0
        x.temp <- x.obs * cos(theta) + y.obs * sin(theta)
        y.temp <- y.obs * cos(theta) - x.obs * sin(theta)
        x.obs  <- x.temp
        y.obs  <- y.temp
        x.theo <- x.obs
        y.theo <- expr(z[4:p], x.obs)
        Value  <- sum((y.obs - y.theo)^2)
        return(Value)
    }
    for (i in 1:nrow(ini.val)) {
      if(method=="Nelder-Mead"){
        res <- optim(ini.val[i, ], obj.fun, method=method, control = control)
      }
      if(method=="L-BFGS-B"){
        res <- optim(ini.val[i, ], obj.fun, method=method, 
                     lower=lower, upper=upper, control = control)
      }
      mat[i, ] <- c(res$par, res$val)
    }
    Names <- rep(NA, len = p)
    for (k in 1:p) {
        Names[k] <- paste("P[", k, "]", sep = "")
    }
    colnames(mat) <- c(Names, "RSS")
    ind <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
    MAT <- mat
    par <- as.vector(mat[ind, 1:p])
    PAR <- par[4:p]
    if (par[3] > 2 * pi) {
        par[3] <- par[3]%%(2 * pi)
    }
    if (par[3] < -2 * pi) {
        par[3] <- -((-par[3])%%(2 * pi))
    }
    goal.x0    <- par[1]
    goal.y0    <- par[2]
    goal.theta <- par[3]
    x.obs   <- x - goal.x0
    y.obs   <- y - goal.y0
    x.new   <- x.obs * cos(goal.theta) + y.obs * sin(goal.theta)
    y.new   <- y.obs * cos(goal.theta) - x.obs * sin(goal.theta)
    x.obs   <- x.new
    y.obs   <- y.new
    x1      <- x.obs
    y1      <- y.obs
    x.theo  <- x.obs
    y.theo  <- expr(PAR, x.obs)
    x2      <- x.theo
    y2      <- y.theo
    r.obs   <- sqrt(x1^2 + y1^2)
    phi.obs <- acos(x1/r.obs)
    phiA    <- phi.obs

    if (is.null(angle)) {
        x.new    <- x1 * cos(goal.theta) - y1 * sin(goal.theta) + par[1]
        y.new    <- y1 * cos(goal.theta) + x1 * sin(goal.theta) + par[2]
        x.theory <- x2 * cos(goal.theta) - y2 * sin(goal.theta) + par[1]
        y.theory <- y2 * cos(goal.theta) + x2 * sin(goal.theta) + par[2]
        phiB     <- phiA + goal.theta
    }
    if (!is.null(angle)) {
        par.new <- par
        par.new[3] <- angle
        x.new    <- x1 * cos(angle) - y1 * sin(angle) + par.new[1]
        y.new    <- y1 * cos(angle) + x1 * sin(angle) + par.new[2]
        x.theory <- x2 * cos(angle) - y2 * sin(angle) + par[1]
        y.theory <- y2 * cos(angle) + x2 * sin(angle) + par[2]
        phiB     <- phiA + angle
    }
    if (is.null(xlim)) 
        xlim <- NULL
    if (is.null(ylim)) 
        ylim <- NULL
    if (!is.null(xlim)) 
        xlim <- xlim
    if (!is.null(ylim)) 
        ylim <- ylim
    if (!is.null(unit)) {
        xlabel <- bquote(paste(italic("x"), " (", .(unit), ")", 
            sep = ""))
        ylabel <- bquote(paste(italic("y"), " (", .(unit), ")", 
            sep = ""))
    }
    if (is.null(unit)) {
        xlabel <- bquote(italic("x"))
        ylabel <- bquote(italic("y"))
    }
    if (stand.fig == "T" | stand.fig == "TRUE" | stand.fig == 
        "True") {
        dev.new()
        par( mar=c(5, 5 ,2, 2) )
        par( mgp=c(3, 1, 0) )  
        plot(c(x1, x2), c(y1, y2), xlab = xlabel, ylab = ylabel, las = 1, cex.lab = 1.5, 
            cex.axis = 1.5, type = "n", asp = 1)
        points(x1, y1, cex = 1, col = "grey40")
        lines(x2, y2, col = 2, lwd = 2)
        title(main = main, cex.main = 1.5, col.main = 4, font.main = 1)
        abline(h = 0, lty = 2, col = 4)
        abline(v = 0, lty = 2, col = 4)
    }
    if(fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True") {
        x_ran <- c(x.new, x.theory)
        y_ran <- c(y.new, y.theory)
        if (is.null(angle)) {
            dev.new()
            par( mar=c(5, 5 ,2, 2) )
            par( mgp=c(3, 1, 0) )  
            plot(x_ran, y_ran, xlab = xlabel, ylab = ylabel, 
                type = "n", asp = 1, xlim = xlim, ylim = ylim, 
                cex.lab = 1.5, cex.axis = 1.5)
            points(x.new, y.new, cex = 1, col = "grey40")
            lines(x.theory, y.theory, col = 2, lwd = 2)

            title(main = main, cex.main = 1.5, col.main = 4, 
                font.main = 1)
            if (abs(par[3])%%pi/2 != 0) {
                slope <- tan(par[3])
                abline(-slope * par[1] + par[2], slope, col = 1, lty = 2)
                abline(v = par[1], col = 4, lty = 2, lwd = 1)
                abline(h = par[2], col = 4, lty = 2, lwd = 1)
            }
            if (abs(par[3])%%pi/2 == 0) {
                abline(v = par[1], col = 1, lty = 2)
                abline(v = par[1], col = 4, lty = 2, lwd = 1)
                abline(h = par[2], col = 4, lty = 2, lwd = 1)
            }
        }
        if (!is.null(angle)) {
            dev.new()
            par( mar=c(5, 5 ,2, 2) )
            par( mgp=c(3, 1, 0) )  
            plot(x_ran, y_ran, xlab = xlabel, ylab = ylabel, 
                type = "n", asp = 1, xlim = xlim, ylim = ylim, 
                cex.lab = 1.5, cex.axis = 1.5)
            points(x.new, y.new, cex = 1, col = "grey40")
            lines(x.theory, y.theory, col = 2, lwd = 2)

            title(main = main, cex.main = 1.5, col.main = 4, 
                font.main = 1)
            if (abs(par.new[3])%%pi/2 != 0) {
                slope <- tan(par.new[3])
                abline(-slope * par[1] + par[2], slope, col = 1, lty = 2)
                abline(v = par[1], col = 4, lty = 2, lwd = 1)
                abline(h = par[2], col = 4, lty = 2, lwd = 1)
            }
            if (abs(par.new[3])%%pi/2 == 0) {
                abline(v = par[1], col = 1, lty = 2)
                abline(v = par[1], col = 4, lty = 2, lwd = 1)
                abline(h = par[2], col = 4, lty = 2, lwd = 1)
            }
        }
    }

    RSS  <- sum((y.obs - y.theo)^2)
    r.sq <- 1 - sum((y.obs - y.theo)^2)/sum((y.obs - mean(y.obs))^2)

    para.tab <- data.frame(Parameter = c(Names, "r.sq", "RSS", 
        "sample.size"), Estimate = c(par, r.sq, RSS, length(x)))
    if(par.list == "T" | par.list == "TRUE" | par.list == "True"){
        print(para.tab)
        cat("\n")
    }
    return(list(par = par, r.sq = r.sq, RSS = RSS, sample.size = length(x), 
           x.stand.obs = x1, x.stand.pred = x2, y.stand.obs = y1, 
           y.stand.pred = y2, x.obs = x.new, x.pred = x.theory, 
           y.obs = y.new, y.pred = y.theory))
}
