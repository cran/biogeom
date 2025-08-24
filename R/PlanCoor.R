PlanCoor <- function( folder.name, lower.val = 0, upper.val = 250, ratio = 1,  
                fig.opt = TRUE, np = NULL, unit = "unitless", verbose = TRUE ){

  if(!is.null(np) & !is.numeric(np)) 
    stop("'np' should be a integer or NULL!")
  if(!is.null(np) & is.numeric(np)){ 
    if(np <= 0)
      stop("'np' should be a positive integer!")
  }

  if(!is.character(folder.name))
    stop("'folder.name' should be an object of type 'character'!")

  files    <- list.files(folder.name) 
  U        <- length(files)
  FileName <- c()
  x_vec    <- c()
  y_vec    <- c()

  for(counter in 1:U){
    if(verbose){
      print(paste("The current progress is ", counter, "/", U, sep=""))
    }
    x.fig  <- files[counter]
    y.fig  <- substr(x.fig, 1, nchar(x.fig)-4)
    Format <- substr(x.fig, nchar(x.fig)-2, nchar(x.fig))
    Format <- tolower(Format)
    if( Format != "bmp" )
      stop("All files in this folder should be grayscale .bmp files!")
    if( Format == "bmp" ){
      suppressWarnings({
        mat <- read.bmp( paste(folder.name, "/", files[counter], sep="") )
      })
    }

    #### To convert arrays of pixel data between different display conventions ####
    smat <- transmat(mat, from="European", to="spatstat")
    ###############################################################################

    temp.w   <- ncol(mat)
    temp.h   <- nrow(mat)
    temp.res <- attr(mat, "header")$hres/100
    Width    <- temp.w/temp.res * ratio
    Height   <- temp.h/temp.res * ratio

    vim <- im( smat, xrange=c(0, Width), yrange=c(0, Height),
               unitname=unit )
    xlabel <- bquote( paste(italic(x), " (", .(unit), ")", sep="") ) 
    ylabel <- bquote( paste(italic(y), " (", .(unit), ")", sep="") )

    X   <- data.frame(vim)
    x1  <- X$x
    y1  <- X$y
    z1  <- X$value  
    ind <- z1 >= lower.val & z1 <= upper.val

    x   <- x1[ind]
    y   <- y1[ind]
    z   <- z1[ind]

    if(!is.null(np)){
      if(length(x) > np){
          ind2 <- sample(1:length(x), np, replace=FALSE)
          x <- x[ind2]
          y <- y[ind2]
      }
    }

    if( fig.opt == "TRUE" ){
      dev.new()
      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )
      par(family="serif")
      plot( vim, cex.lab=1.5, cex.axis=1.5, main="", axes=TRUE )
      mtext(xlabel, side = 1, line = 2, cex=1.5)
      mtext(ylabel, side = 2, line = 2, cex=1.5)
      title(main=y.fig, cex.main=1.5, col.main=4)

      dev.new()
      par( mar=c(5, 5, 2, 2) )
      par( mgp=c(3.5, 1, 0) )
      par(family="serif")
      plot( x, y, asp=1, cex.lab=1.5, cex.axis=1.5, main="", 
            xlab=xlabel, ylab=ylabel, xlim=c(0, Width), ylim=c(0, Height) )
      title(main=y.fig, cex.main=1.5, col.main=4)
    }
  
    x_vec    <- c(x_vec, x)
    y_vec    <- c(y_vec, y)
    FileName <- c(FileName, rep(y.fig, len=length(x))) 

  }

  return( list(FileName=FileName, x=x_vec, y=y_vec) )
}







