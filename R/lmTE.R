lmTE <- function (x, y, dev.angle = NULL,  
    weights = NULL, fig.opt = TRUE, prog.opt = TRUE, xlim = NULL, 
    ylim = NULL, unit = NULL, main = NULL) 
{
    w <- weights
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    if (length(x) != length(y)) 
        stop("'x' should have the same data length as 'y'!")
    if(!is.null(dev.angle) & !is.numeric(dev.angle)){
        stop("'dev.angle' should be NULL or a numerical value or a numerical vector!")
    }

    Tem <- cbind(x, y)
    Tem <- na.omit(Tem)
    x   <- Tem[, 1]
    y   <- Tem[, 2]
    poly1 <- owin(xrange = c(min(x)[1], max(x)[1]), 
                  yrange = c(min(y)[1], max(y)[1]))
    poly0 <- as.polygonal(owin(poly = list(x = x, y = y)))
    data0 <- ppp(x, y, window = poly1)
    area0 <- area.owin(poly0)
    perimeter0 <- perimeter(poly0)
    mat <- pairdist(data0)
    qq  <- cbind(which(mat == max(mat)[1], 
               arr.ind = TRUE)[1,])
    q <- as.numeric(qq)
    x1 <- x[q[1]]
    y1 <- y[q[1]]
    x2 <- x[q[2]]
    y2 <- y[q[2]]
    xx <- x - x2
    yy <- y - y2
    x1 <- x1 - x2
    y1 <- y1 - y2
    x2 <- x2 - x2
    y2 <- y2 - y2
    z  <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    theta <- acos((x1 - x2)/z)%%(2 * pi)
    if(theta > 2*pi) theta <- theta %% (2 * pi)
    xx2     <- xx * cos(theta) + yy * sin(theta)
    yy2     <- yy * cos(theta) - xx * sin(theta)
    length1 <- max(xx2) - min(xx2)
    xx3     <- xx2 - min(xx2)
    yy3     <- yy2 - min(yy2)
    Index1  <- which.max(yy3)
    Index2  <- which.min(yy3)
    x.mb    <- (xx3[Index1]+xx3[Index2])/2
    if(x.mb > length1/2){ 
      theta <- (theta + pi) %% (2 * pi)
    }
    if (is.null(dev.angle)) {    
      x.new <- xx * cos(theta) + yy * sin(theta)
      y.new <- yy * cos(theta) - xx * sin(theta)
      length0 <- max(x.new) - min(x.new)
      width0  <- max(y.new) - min(y.new)
      Index1 <- which.max(y.new)
      Index2 <- which.min(y.new)
      x.mb   <- (x.new[Index1]+x.new[Index2])/2
      if(x.mb > length0/2){         
        xT    <- x.new * cos(pi) + y.new * sin(pi)
        yT    <- y.new * cos(pi) - x.new * sin(pi)
        x.new <- xT
        y.new <- yT
        theta <- (theta + pi) %% (2 * pi)
        x.new <- x.new - min(x.new)
      }
      x.new  <- x.new - (max(x.new)+min(x.new))/2      
      y.new  <- y.new - (max(y.new)+min(y.new))/2
      xv     <- x.new/(length0/2)
      yv     <- y.new/(length0/2)
      xv     <<- xv
      yv     <<- yv
      inde1  <- which(yv > 0)
      inde2  <- which(yv <= 0)
      log.yv <- yv
      log.yv[inde1] <- log(yv[inde1])
      log.yv[inde2] <- log(-yv[inde2])
      temp <- 1 - xv^2
      ind0 <- which(temp <= 0)
      temp[ind0] <- NA
      z1   <- xv
      z2   <- xv^2
      z3   <- log(temp)
      res  <- lm(log.yv ~ z1 + z2 + offset(0.5*z3), weights = w) 
      res  <<- res
      alpha0 <- as.numeric( res$coefficients[[1]] )
      alpha1 <- as.numeric( res$coefficients[[2]] )
      alpha2 <- as.numeric( res$coefficients[[3]] )
      yv.pred <- TE(P = c(alpha0, alpha1, alpha2 ), x = xv)
      yv.pred[inde2] <- -yv.pred[inde2]
      y.pred <- yv.pred * length0/2 
      RSS1   <- sum((yv - yv.pred)^2)
      RSS2   <- sum((y.new - y.pred)^2)
      RMSE1  <- sqrt(sum((yv - yv.pred)^2)/length(yv))
      RMSE2  <- sqrt(sum((y.new - y.pred)^2)/length(y.new))
    }

    if (!is.null(dev.angle)) {   
      RSSV <- c()
      AngleTotal <- length(dev.angle) + 1
      for(i in 1:length(dev.angle)){
        if(prog.opt == "TRUE" | prog.opt == "T" | prog.opt == "True")
            print(paste("Progress: ", i, "/", AngleTotal, sep=""))
        epsilon <- dev.angle[i]
        x.new   <- xx * cos(theta+epsilon) + yy * sin(theta+epsilon)
        y.new   <- yy * cos(theta+epsilon) - xx * sin(theta+epsilon)
        length0 <- max(x.new) - min(x.new)
        width0  <- max(y.new) - min(y.new)  
        Index1 <- which.max(y.new)
        Index2 <- which.min(y.new)
        x.mb   <- (x.new[Index1]+x.new[Index2])/2
        if(x.mb > length0/2){         
          xT    <- x.new * cos(pi) + y.new * sin(pi)
          yT    <- y.new * cos(pi) - x.new * sin(pi)
          x.new <- xT
          y.new <- yT
          theta <- (theta + pi) %% (2 * pi)
          x.new <- x.new - min(x.new)
        }
        x.new  <- x.new - (max(x.new)+min(x.new))/2      
        y.new  <- y.new - (max(y.new)+min(y.new))/2
        xv     <- x.new/(length0/2)
        yv     <- y.new/(length0/2)
        inde1  <- which(yv > 0)
        inde2  <- which(yv <= 0)
        log.yv <- yv
        log.yv[inde1] <- log(yv[inde1])
        log.yv[inde2] <- log(-yv[inde2])
        temp <- 1 - xv^2
        ind0 <- which(temp <= 0)
        temp[ind0] <- NA
        z1   <- xv
        z2   <- xv^2
        z3   <- log(temp)
        z1[inde2] <- xv[inde2]
        z2[inde2] <- xv[inde2]^2
        z3[inde2] <- log(temp[inde2])  
        res  <- lm(log.yv ~ z1 + z2 + offset(0.5*z3), weights = w) 
        res  <<- res
        alpha0 <- as.numeric( res$coefficients[[1]] )
        alpha1 <- as.numeric( res$coefficients[[2]] )
        alpha2 <- as.numeric( res$coefficients[[3]] )
        yv.pred <- TE(P = c(alpha0, alpha1, alpha2), x = xv)
        yv.pred[inde2] <- -yv.pred[inde2]
        y.pred <- yv.pred * length0/2 
        RSS1   <- sum((yv - yv.pred)^2)
        RSSV   <- c(RSSV, RSS1)  
      }
      if(prog.opt == "TRUE" | prog.opt == "T" | prog.opt == "True")
          print(paste("Progress: ", AngleTotal, "/", AngleTotal, sep=""))
      Ind0     <- which.min(RSSV)[1]
      epsilon0 <- dev.angle[Ind0]
      x.new <- xx * cos(theta+epsilon0) + yy * sin(theta+epsilon0)
      y.new <- yy * cos(theta+epsilon0) - xx * sin(theta+epsilon0)
      length0 <- max(x.new) - min(x.new)
      width0  <- max(y.new) - min(y.new)
      Index1 <- which.max(y.new)
      Index2 <- which.min(y.new)
      x.mb   <- (x.new[Index1]+x.new[Index2])/2
      if(x.mb > length0/2){         
        xT    <- x.new * cos(pi) + y.new * sin(pi)
        yT    <- y.new * cos(pi) - x.new * sin(pi)
        x.new <- xT
        y.new <- yT
        theta <- (theta + pi) %% (2 * pi)
        x.new <- x.new - min(x.new)
      }
      x.new  <- x.new - (max(x.new)+min(x.new))/2      
      y.new  <- y.new - (max(y.new)+min(y.new))/2
      xv     <- x.new/(length0/2)
      yv     <- y.new/(length0/2)
      inde1  <- which(yv > 0)
      inde2  <- which(yv <= 0)
      log.yv <- yv
      log.yv[inde1] <- log(yv[inde1])
      log.yv[inde2] <- log(-yv[inde2])
      temp <- 1 - xv^2
      ind0 <- which(temp <= 0)
      temp[ind0] <- NA
      z1   <- xv
      z2   <- xv^2
      z3   <- log(temp)
      z1[inde2] <- xv[inde2]
      z2[inde2] <- xv[inde2]^2
      z3[inde2] <- log(temp[inde2])  
      res  <- lm(log.yv ~ z1 + z2 + offset(0.5*z3), weights = w) 
      res  <<- res
      alpha0 <- as.numeric( res$coefficients[[1]] )
      alpha1 <- as.numeric( res$coefficients[[2]] )
      alpha2 <- as.numeric( res$coefficients[[3]] )
      yv.pred <- TE(P = c(alpha0, alpha1, alpha2), x = xv)
      yv.pred[inde2] <- -yv.pred[inde2]
      y.pred <- yv.pred * length0/2 
      RSS1   <- sum((yv - yv.pred)^2)
      RSS2   <- sum((y.new - y.pred)^2)
      RMSE1  <- sqrt(sum((yv - yv.pred)^2)/length(yv))
      RMSE2  <- sqrt(sum((y.new - y.pred)^2)/length(y.new))
    }

    if (!is.null(unit)) {
        xlabel <- bquote(paste(italic(x), " (", .(unit), ")", 
            sep = ""))
        ylabel <- bquote(paste(italic(y), " (", .(unit), ")", 
            sep = ""))
    }
    if (is.null(unit)) {
        xlabel <- bquote(italic(x))
        ylabel <- bquote(italic(y))
    }
    if (is.null(xlim)) 
        xlim <- NULL
    if (is.null(ylim)) 
        ylim <- NULL
    if (!is.null(xlim)) 
        xlim <- xlim
    if (!is.null(ylim)) 
        ylim <- ylim
    if (fig.opt == "TRUE" | fig.opt == "T" | fig.opt == "True") {
        dev.new()
        par(family = "serif")
        par(mar = c(5, 5, 2, 2))
        plot(xv, yv, asp = 1, col = "grey50", lwd = 2, xlab = bquote(italic(x)), 
            ylab = bquote(italic(y)), xlim = xlim, ylim = ylim, 
            cex.lab = 1.5, cex.axis = 1.5, type = "l")
        lines(xv, yv.pred, col = 2)
        title(main = main, cex.main = 1.5, col.main = 4, font.main = 1)
        text(0, 0, bquote(paste("RMSE = ", .(round(RMSE1, 4)), 
            sep = "")), pos = NULL, cex = 1.5, col = 1)
        dev.new()
        par(family = "serif")
        par(mar = c(5, 5, 2, 2))
        plot(x.new, y.new, asp = 1, col = "grey50", lwd = 2, 
            xlab = xlabel, ylab = ylabel, xlim = xlim, ylim = ylim, 
            cex.lab = 1.5, cex.axis = 1.5, type = "l")
        lines(x.new, y.pred, col = 2)
        title(main = main, cex.main = 1.5, col.main = 4, font.main = 1)
        text(0, 0, bquote(paste("RMSE = ", .(round(RMSE2, 4)), 
            sep = "")), pos = NULL, cex = 1.5, col = 1)
    }

    if (is.null(dev.angle)){
      Res <- list(lm.te = res, par = c(alpha0, alpha1, alpha2), theta = theta,          
          x.obs = x.new, y.obs = y.new, y.pred = y.pred, x.stand.obs = xv, 
          y.stand.obs = yv, y.stand.pred = yv.pred, scan.length = length0, 
          scan.width = width0, scan.area = area0, scan.perimeter = perimeter0, 
          RSS.scaled = RSS1, RSS = RSS2, sample.size = length(y.new), 
          RMSE.scaled = RMSE1, RMSE = RMSE2)
    }

    if (!is.null(dev.angle)){
      Res <- list(lm.te = res, par = c(alpha0, alpha1, alpha2), theta = theta, 
          epsilon = epsilon0, RSS.vector = RSSV, 
          x.obs = x.new, y.obs = y.new, y.pred = y.pred, x.stand.obs = xv, 
          y.stand.obs = yv, y.stand.pred = yv.pred, scan.length = length0, 
          scan.width = width0, scan.area = area0, scan.perimeter = perimeter0, 
          RSS.scaled = RSS1, RSS = RSS2, sample.size = length(y.new), 
          RMSE.scaled = RMSE1, RMSE = RMSE2)
    }
    return(Res)
}
