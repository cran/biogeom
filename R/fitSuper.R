
fitSuper <- function (x, y, ini.val,  
    control = list(), par.list = FALSE, stand.fig = TRUE, angle = NULL, 
    fig.opt = FALSE, np = 2000, xlim = NULL, ylim = NULL, unit = NULL, 
    main = NULL){
    if (length(x) != length(y)) 
        stop("'x' should have the same data length as 'y'!")
    Tem <- cbind(x, y)
    Tem <- na.omit(Tem)
    x <- Tem[, 1]
    y <- Tem[, 2]
    ini.val <- as.list(ini.val)
    p <- length(ini.val)
    s <- 1
    for (i in 1:p) {
        s <- s * length(ini.val[[i]])
    }
    ini.val <- expand.grid(ini.val)
    mat <- matrix(NA, nrow = s, ncol = (p + 1))
    x.obs <- x
    y.obs <- y
    obj.fun <- function(z) {
        x0 <- z[1]
        y0 <- z[2]
        theta <- z[3]
        x.obs <- x.obs - x0
        y.obs <- y.obs - y0
        x.temp <- x.obs * cos(theta) + y.obs * sin(theta)
        y.temp <- y.obs * cos(theta) - x.obs * sin(theta)
        x.obs <- x.temp
        y.obs <- y.temp
        r.obs <- sqrt(x.obs^2 + y.obs^2)
        phi.obs <- acos(x.obs/r.obs)
        cond1 <- y.obs >= 0
        cond2 <- y.obs < 0
        xx1 <- x.obs[cond1]
        xx2 <- x.obs[cond2]
        yy1 <- y.obs[cond1]
        yy2 <- y.obs[cond2]
        phi1 <- phi.obs[cond1]
        phi2 <- 2 * pi - phi.obs[cond2]
        phi.obs <- c(phi1, phi2)
        x.obs <- c(xx1, xx2)
        y.obs <- c(yy1, yy2)
        resu <- sort(phi.obs, decreasing = FALSE, index.return = T)
        phi.obs <- resu$x
        Index <- resu$ix
        x.obs <- x.obs[Index]
        y.obs <- y.obs[Index]
        r.obs <- sqrt(x.obs^2 + y.obs^2)
        r.theo <- GE(P = z[4:p], phi = phi.obs, m = 4, simpver = 9)
        x.theo <- r.theo * cos(phi.obs)
        y.theo <- r.theo * sin(phi.obs)
        temp   <- sum((r.obs - r.theo)^2)
        if(z[5] > 1 | z[5] < 0) temp <- Inf    
        return(temp)
    }
    for (i in 1:nrow(ini.val)) {
        res <- optim(ini.val[i, ], obj.fun, control = control)
        mat[i, ] <- c(res$par, res$val)
    }
    Names <- rep(NA, len = p)
    for (k in 1:p) {
        Names[k] <- paste("P[", k, "]", sep = "")
    }
    colnames(mat) <- c(Names, "RSS")
    ind <- which(mat[, p + 1] == min(mat[, p + 1])[1])[1]
    par <- as.vector(mat[ind, 1:p])
    if (par[3] > 2 * pi) {
        par[3] <- par[3]%%(2 * pi)
    }
    if (par[3] < -2 * pi) {
        par[3] <- -((-par[3])%%(2 * pi))
    }
    goal.x0 <- par[1]
    goal.y0 <- par[2]
    goal.theta <- par[3]
    x.obs <- x - goal.x0
    y.obs <- y - goal.y0
    x.new <- x.obs * cos(goal.theta) + y.obs * sin(goal.theta)
    y.new <- y.obs * cos(goal.theta) - x.obs * sin(goal.theta)
    x.obs <- x.new
    y.obs <- y.new
    r.obs <- sqrt(x.obs^2 + y.obs^2)
    phi.obs <- acos(x.obs/r.obs)
    cond1 <- y.obs >= 0
    cond2 <- y.obs < 0
    xx1  <- x.obs[cond1]
    xx2  <- x.obs[cond2]
    yy1  <- y.obs[cond1]
    yy2  <- y.obs[cond2]
    phi1 <- phi.obs[cond1]
    phi2 <- 2 * pi - phi.obs[cond2]
    phi.obs <- c(phi1, phi2)
    x.obs <- c(xx1, xx2)
    y.obs <- c(yy1, yy2)
    resu  <- sort(phi.obs, decreasing = FALSE, index.return = T)
    phi.obs <- resu$x
    Index <- resu$ix
    x.obs <- x.obs[Index]
    y.obs <- y.obs[Index]
    r.obs <- sqrt(x.obs^2 + y.obs^2)
    x1 <- x.obs
    y1 <- y.obs
    poly0 <- as.polygonal(owin(poly = list(x = x1, y = y1)))
    Area  <- area(poly0)
    Perimeter <- perimeter(poly0)
    r.theo <- GE(par[4:p], phi = phi.obs, m = 4, simpver = 9)
    x.theo <- r.theo * cos(phi.obs)
    y.theo <- r.theo * sin(phi.obs)
    x2 <- x.theo
    y2 <- y.theo
    phiA <- phi.obs
    phi.arr <- seq(0, 2 * pi, len = np)
    if (is.null(angle)) {
        results <- curveGE(GE, P = par, phi = phi.obs, m = 4, 
                           simpver = 9)
        results2 <- curveGE(GE, P = par, phi = phi.arr, m = 4, 
                            simpver = 9)
        x.new <- x1 * cos(goal.theta) - y1 * sin(goal.theta) + 
            par[1]
        y.new <- y1 * cos(goal.theta) + x1 * sin(goal.theta) + 
            par[2]
        r.new <- sqrt((x.new - par[1])^2 + (y.new - par[2])^2)
        phiB <- phiA + goal.theta
    }
    if (!is.null(angle)) {
        par.new <- par
        par.new[3] <- angle
        results <- curveGE(GE, P = par.new, phi = phi.obs, 
                           m = 4, simpver = 9)
        results2 <- curveGE(GE, P = par.new, phi = phi.arr, 
                            m = 4, simpver = 9)
        x.new <- x1 * cos(angle) - y1 * sin(angle) + par.new[1]
        y.new <- y1 * cos(angle) + x1 * sin(angle) + par.new[2]
        r.new <- sqrt((x.new - par.new[1])^2 + (y.new - par.new[2])^2)
        phiB <- phiA + angle
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
        plot(x1, y1, xlab = xlabel, ylab = ylabel, cex.lab = 1.5, 
            cex.axis = 1.5, type = "l", lwd = 3, col = "grey50", 
            asp = 1)
        lines(x2, y2, col = 2, lwd = 2)
        title(main = main, cex.main = 1.5, col.main = 4, font.main = 1)
        abline(h = 0, lty = 2, col = 4)
        abline(v = 0, lty = 2, col = 4)
    }
    if (fig.opt == "T" | fig.opt == "TRUE" | fig.opt == "True") {
        if (is.null(angle)) {
            dev.new()
            plot(x.new, y.new, xlab = xlabel, ylab = ylabel, 
                type = "l", cex.lab = 1.5, cex.axis = 1.5, col = "grey50", 
                lwd = 3, asp = 1, xlim = xlim, ylim = ylim)
            lines(results2$x, results2$y, type = "l", asp = 1, 
                col = 2, lwd = 2)
            title(main = main, cex.main = 1.5, col.main = 4, 
                font.main = 1)
            if (abs(par[3])%%pi/2 != 0) {
                slope <- tan(par[3])
                abline(-slope * par[1] + par[2], slope, col = 1, 
                  lty = 2)
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
            plot(x.new, y.new, xlab = xlabel, ylab = ylabel, 
                type = "l", cex.lab = 1.5, cex.axis = 1.5, col = "grey50", 
                lwd = 3, asp = 1, xlim = xlim, ylim = ylim)
            lines(results2$x, results2$y, type = "l", asp = 1, 
                col = 2, lwd = 2)
            title(main = main, cex.main = 1.5, col.main = 4, 
                font.main = 1)
            if (abs(par.new[3])%%pi/2 != 0) {
                slope <- tan(par.new[3])
                abline(-slope * par[1] + par[2], slope, col = 1, 
                  lty = 2)
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
    x.range <- range(x1)
    y.range <- range(y1)
    Length  <- x.range[2] - x.range[1]
    Width   <- y.range[2] - y.range[1]
    r1      <- sqrt(x1^2 + y1^2)
    r2      <- sqrt(x2^2 + y2^2)
    RSS     <- sum((r1 - r2)^2)
    r.sq <- 1 - sum((r1 - r2)^2)/sum((r1 - mean(r1))^2)
    para.tab <- data.frame(Parameter = c(Names, "Length", "Width", 
        "r.sq", "RSS", "sample.size"), Estimate = c(par, Length, 
        Width, r.sq, RSS, length(x)))
    if (par.list == "T" | par.list == "TRUE" | par.list == "True") {
        print(para.tab)
        cat("\n")
    }
    return(list(par = par, scan.length = Length, scan.width = Width, 
        scan.area = Area, scan.perimeter = Perimeter, r.sq = r.sq, 
        RSS = RSS, sample.size = length(x), 
        phi.stand.obs = phiA, phi.trans = phiB, r.stand.obs = r1, 
        r.stand.pred = r2, x.stand.obs = x1, x.stand.pred = x2, 
        y.stand.obs = y1, y.stand.pred = y2, r.obs = r.new, r.pred = results$r, 
        x.obs = x.new, x.pred = results$x, y.obs = y.new, y.pred = results$y))
}

