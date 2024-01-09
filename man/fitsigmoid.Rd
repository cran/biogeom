\name{fitsigmoid}
\encoding{UTF-8}
\alias{fitsigmoid}
\title{
Data-Fitting Function for the Sigmoid Growth Equation
}
\description{
\code{fitsigmoid} is used to estimate the parameters of a sigmoid growth equation based on the integral of 
a performance equation or one of its simplified versions.
}
\usage{
fitsigmoid(expr, x, y, ini.val, simpver = 1, 
           control = list(), par.list = FALSE, fig.opt = FALSE,
           xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,         
           main = NULL, subdivisions = 100L,
           rel.tol = .Machine$double.eps^0.25, 
           abs.tol = rel.tol, stop.on.error = TRUE, 
           keep.xy = FALSE, aux = NULL)
}
\arguments{
  \item{expr}{a performance equation or one of its simplified versions that is used to build a sigmoid growth equation.}
  \item{x}{the observed investigation times.}
  \item{y}{the observed \eqn{y} values (i.e., biomass, height, body length, etc.).}
  \item{ini.val}{the initial values of the model parameters.}
  \item{simpver}{an optional argument to use the simplified version of the performance equation.}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{par.list}{the option of showing the list of parameters on the screen.}
  \item{fig.opt}{an optional argument to draw the observations and the predicted sigmoid curve.}
  \item{xlim}{the range of the \eqn{x}-axis over which to plot a sigmoid growth curve.}
  \item{ylim}{the range of the \eqn{y}-axis over which to plot a sigmoid growth curve.}
  \item{xlab}{the label of the \eqn{x}-axis when showing a sigmoid growth curve.}
  \item{ylab}{the label of the \eqn{y}-axis when showing a sigmoid growth curve.}
  \item{main}{the main title of the figure.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
}
\details{
Here, \code{ini.val} only includes the initial values of the model parameters as a list. 
The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization of minimizing the residual 
sum of squares (RSS) between the observed and predicted \eqn{y} values. The \code{\link{optim}} 
function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
The performance equations denote \code{\link{MbetaE}}, \code{\link{MBriereE}}, 
\code{\link{MLRFE}}, \code{\link{MPerformanceE}} and their simplified versions. 
The arguments of \code{P} and \code{simpver} should correspond 
to \code{expr} (i.e., \code{\link{MbetaE}} or \code{\link{MBriereE}} 
or \code{\link{MLRFE}} or \code{\link{MPerformanceE}}). 
The sigmoid equation is the integral of a performance equation or one of its simplified versions.
}
\value{
  \item{par}{the estimates of the model parameters.}
  \item{r.sq}{the coefficient of determination between the observed and predicted \eqn{y} values.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values.}
  \item{sample.size}{the number of data points used in the data fitting.}
  \item{x}{the observed \eqn{x} values.}
  \item{y}{the observed \eqn{y} values.}
  \item{y.pred}{the predicted \eqn{y} values.}       
}
\note{
Here, the user can define other performance equations, but new equations or 
their simplified versions should include the lower and upper thresholds on 
the \eqn{x}-axis corresponding to \eqn{y = 0}, whose indices should 
be the same as those in \code{\link{MbetaE}} or \code{\link{MBriereE}} 
or \code{\link{MLRFE}} or \code{\link{MPerformanceE}}.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Jin, J., Quinn, B.K., Shi, P. (2022) The modified Brière equation and its 
applications. \emph{Plants} 11, 1769. \doi{10.3390/plants11131769}

Lian, M., Shi, P., Zhang, L., Yao, W., Gielis, J., Niklas, K.J. (2023) A generalized performance equation 
and its application in measuring the Gini index of leaf size inequality. 
\emph{Trees \eqn{-} Structure and Function} 37, 1555\eqn{-}1565. \doi{10.1007/s00468-023-02448-8}

Nelder, J.A., Mead, R. (1965) A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Fan, M., Ratkowsky, D.A., Huang, J., Wu, H., Chen, L., Fang, S., 
Zhang, C. (2017) Comparison of two ontogenetic growth equations for animals and plants. 
\emph{Ecological Modelling} 349, 1\eqn{-}10. \doi{10.1016/j.ecolmodel.2017.01.012}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}
}
\seealso{
\code{\link{areaovate}}, \code{\link{MbetaE}}, \code{\link{MBriereE}}, \code{\link{MLRFE}}, 
\code{\link{MPerformanceE}}, \code{\link{sigmoid}}
}
\examples{

# The shrimp growth data(See the supplementary table in West et al., 2001)
# West, G.B., Brown, J.H., Enquist, B.J. (2001) A general model for ontogenetic growth. 
#     Nature 413, 628-631.
t0 <- c(3, 60, 90, 120, 150, 180, 384)
m0 <- c(0.001, 0.005, 0.018, 0.037, 0.06, 0.067, 0.07)

dev.new()
plot( t0, m0, cex.lab=1.5, cex.axis=1.5, col=4,
      xlab=expression(italic(x)), ylab=expression(italic(y)) )

xopt0   <- seq(100, 150, by=5)
ini.val <- list(0.035, xopt0, 200, 1)
resu1   <- fitsigmoid(MLRFE, x=t0, y=m0, ini.val=ini.val, simpver=1, fig.opt=TRUE, par.list=TRUE)

\donttest{
  xopt0   <- seq(100, 150, by=5)
  ini.val <- list(0.035, xopt0, 200, 1)
  resu1   <- fitsigmoid(MbetaE, x=t0, y=m0, ini.val=ini.val, simpver=1, fig.opt=TRUE)

  m.ini   <- c(0.5, 1, 2, 3, 4, 5, 10, 20)
  ini.val <- list(1e-8, m.ini, 200, 1)
  resu2   <- fitsigmoid(MBriereE, x=t0, y=m0, ini.val=ini.val, simpver=1, 
             fig.opt=TRUE, control=list(reltol=1e-20, maxit=20000, trace=FALSE), 
             subdivisions=100L, rel.tol=.Machine$double.eps^0.25, 
             abs.tol=.Machine$double.eps^0.25, stop.on.error=TRUE, 
             keep.xy=FALSE, aux=NULL)

  delta0  <- c(0.5, 1, 2, 5, 10, 20)
  ini.val <- list(0.035, 150, -100, 200, delta0)
  resu3   <- fitsigmoid(MLRFE, x=t0, y=m0, ini.val=ini.val, simpver=NULL, 
               fig.opt=TRUE, control=list(reltol=1e-20, maxit=20000), 
               subdivisions = 100L, rel.tol=.Machine$double.eps^0.25, 
               abs.tol=.Machine$double.eps^0.25, stop.on.error=TRUE, 
               keep.xy=FALSE, aux=NULL)

  a.ini   <- c(0.1, 1, 10, 100, 200)
  b.ini   <- 200
  ini.val <- list(0.001, 0.02, 0.15, 0, 200, a.ini, b.ini)
  resu4   <- fitsigmoid(MPerformanceE, x=t0, y=m0, ini.val=ini.val, simpver=NULL, 
             fig.opt=TRUE, control=list(reltol=1e-20, maxit=20000, trace=FALSE), 
             subdivisions=100L, rel.tol=.Machine$double.eps^0.25, 
             abs.tol=.Machine$double.eps^0.25, stop.on.error=TRUE, 
             keep.xy=FALSE, aux=NULL)
  resu5   <- fitsigmoid(MPerformanceE, x=t0, y=m0, ini.val=resu4$par, simpver=NULL, 
               fig.opt=TRUE, control=list(reltol=1e-30, maxit=200000, trace=FALSE))

  ini.val <- list(0.001, 0.01, c(0.1, 1, 10), 0, 200)
  resu6   <- fitsigmoid(MPerformanceE, x=t0, y=m0, ini.val=ini.val, simpver=2, 
             fig.opt=TRUE, control=list(reltol=1e-20, maxit=20000, trace=FALSE), 
             subdivisions=100L, rel.tol=.Machine$double.eps^0.25, 
             abs.tol=.Machine$double.eps^0.25, stop.on.error=TRUE, 
             keep.xy=FALSE, aux=NULL)
}

graphics.off()
}