\name{fitLorenz}
\encoding{UTF-8}
\alias{fitLorenz}
\title{
Data-Fitting Function for the Rotated and Right-Shifted Lorenz Curve
}
\description{
\code{fitLorenz} is used to estimate the parameters of the rotated and right-shifted Lorenz curve 
using version 4 or 5 of \code{\link{MPerformanceE}}, or the Lorenz equations including 
\code{\link{SarabiaE}}, \code{\link{SCSE}}, and \code{\link{SHE}}.
}
\usage{
fitLorenz(expr, z, ini.val, simpver = 4, 
          control = list(), par.list = FALSE, 
          fig.opt = FALSE, np = 2000, 
          xlab=NULL, ylab=NULL, main = NULL, subdivisions = 100L,
          rel.tol = .Machine$double.eps^0.25, 
          abs.tol = rel.tol, stop.on.error = TRUE, 
          keep.xy = FALSE, aux = NULL, par.limit = TRUE)
}
\arguments{
  \item{expr}{version 4 or 5 of \code{\link{MPerformanceE}}, or the Lorenz equations including 
\code{\link{SarabiaE}}, \code{\link{SCSE}}, and \code{\link{SHE}}.}
  \item{z}{the observations of size distribution (i.e., the household income distribution, the leaf size distribution).}
  \item{ini.val}{the initial values of the model parameters.}
  \item{simpver}{an optional argument to use version 4 or 5 of \code{\link{MPerformanceE}}.}
  \item{control}{the list of control parameters for using the \code{\link{optim}} function in package \pkg{stats}.}
  \item{par.list}{the option of showing the list of parameters on the screen.}
  \item{fig.opt}{an optional argument to draw the original and rotated Lorenz curves.}
  \item{np}{the number of data points to draw the predicted original and rotated Lorenz curves.}
  \item{xlab}{the label of the \eqn{x}-axis when showing the original Lorenz curve.}
  \item{ylab}{the label of the \eqn{y}-axis when showing the original Lorenz curve.}
  \item{main}{the main title of the figure.}
  \item{subdivisions}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{rel.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{abs.tol}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{stop.on.error}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{keep.xy}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{aux}{please see the arguments for the \code{\link{integrate}} function in package \pkg{stats}.}
  \item{par.limit}{an optional argument to limit the numerical ranges of model parameters of the three Lorenz 
equations including \code{\link{SarabiaE}}, \code{\link{SCSE}}, and \code{\link{SHE}}.}
}
\details{
Here, \code{ini.val} only includes the initial values of the model parameters as a list. 
The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization of minimizing the residual 
sum of squares (RSS) between the observed and predicted \eqn{y} values. The \code{\link{optim}} 
function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
Here, versions 4 and 5 of \code{\link{MPerformanceE}} and the Lorenz equations including 
\code{\link{SarabiaE}}, \code{\link{SCSE}}, and \code{\link{SHE}} can be used to fit the rotated 
and right-shifted Lorenz curve.

\eqn{\quad} When \code{simpver = 4}, the simplified version 4 of \code{\link{MPerformanceE}} is selected: 
\deqn{\mbox{if } x \in{\left(0, \ \sqrt{2}\right)},}
\deqn{y = c\left(1-e^{-K_{1}x}\right)^{a}\left(1-e^{K_{2}\left(x-\sqrt{2}\right)}\right)^{b};}
\deqn{\mbox{if } x \notin{\left(0, \ \sqrt{2}\right)},}
\deqn{y = 0.}
There are five elements in \code{P}, representing 
the values of \eqn{c}, \eqn{K_{1}}, \eqn{K_{2}}, \eqn{a}, and \eqn{b}, respectively.

\eqn{\quad} When \code{simpver = 5}, the simplified version 5 of \code{\link{MPerformanceE}} is selected:
\deqn{\mbox{if } x \in{\left(0, \ \sqrt{2}\right)},}
\deqn{y = c\left(1-e^{-K_{1}x}\right)\left(1-e^{K_{2}\left(x-\sqrt{2}\right)}\right);}
\deqn{\mbox{if } x \notin{\left(0, \ \sqrt{2}\right)},}
\deqn{y = 0.}
There are three elements in \code{P}, representing 
the values of \eqn{c}, \eqn{K_{1}}, and \eqn{K_{2}}, respectively.

\eqn{\quad} For the Lorenz functions, the user can define any formulae that follow the below form: 
\code{Lorenz.fun <- function(P, x)\{...\}}, where \code{P} is the vector of parameter(s), \code{x} is the preditor that 
ranges between 0 and 1 representing the cumulative proportion of the number of individuals in a statistical 
unit, and \code{Lorenz.fun} is the name of a Lorenz function defined by the user, which also ranges between 0 and 1 
representing the cumulative 
proportion of the income or size in a statistical unit. This package provides three representative Lorenz 
functions: \code{\link{SarabiaE}}, \code{\link{SCSE}}, and \code{\link{SHE}}.  


\eqn{\quad} Here, the Gini coefficient (GC) is calculated as follows when \code{\link{MPerformanceE}} 
is selected:
\deqn{\mbox{GC} = 2\int_{0}^{\sqrt{2}}y\,dx,}
where \eqn{x} and \eqn{y} are the independent and dependent variables 
in version 4 or 5 of \code{\link{MPerformanceE}}, respectively.

\eqn{\quad} However, the Gini coefficient (GC) is calculated as follows when a Lorenz function, e.g., 
\code{\link{SCSE}}, is selected: 
\deqn{\mbox{GC} = 2\int_{0}^{1}y\,dx,}
where \eqn{x} and \eqn{y} are the independent and dependent variables 
in the Lorenz function, respectively. 

\eqn{\quad} For \code{\link{SarabiaE}} and \code{\link{SHE}}, there are explicit formulae for GC 
(Sarabia, 1997; Sitthiyot and Holasut, 2023). 
}
\value{
  \item{x1}{the cumulative proportion of the number of an entity of interest, 
i.e., the number of households of a city, the number of leaves of a plant.}
  \item{y1}{the cumulative proportion of the size of an entity of interest.}
  \item{x}{the \eqn{x} coordinates of the rotated and right-shifted \code{y1} versus \code{x1}.}
  \item{y}{the \eqn{y} coordinates of the rotated and right-shifted \code{y1} versus \code{x1}.}
  \item{par}{the estimates of the model parameters.}
  \item{r.sq}{the coefficient of determination between the observed and predicted \eqn{y} values.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values.}
  \item{sample.size}{the number of data points used in the data fitting.}
  \item{GC}{the calculated Gini coefficient.}       
}
\note{
When \code{\link{MPerformanceE}} is selected, the estimates of the model parameters 
denote those in \code{\link{MPerformanceE}} rather than 
being obtained by directly fitting the \code{y1} versus \code{x1} data; 
when a Lorenz function is selected, the estimates of the model parameters 
denote those in the Lorenz function. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Huey, R.B., Stevenson, R.D. (1979) Integrating thermal physiology and ecology of ectotherms: 
a discussion of approaches. \emph{American Zoologist} 19, 357\eqn{-}366. \doi{10.1093/icb/19.1.357}

Lian, M., Shi, P., Zhang, L., Yao, W., Gielis, J., Niklas, K.J. (2023) A generalized performance equation 
and its application in measuring the Gini index of leaf size inequality. 
\emph{Trees \eqn{-} Structure and Function} 37, 1555\eqn{-}1565. \doi{10.1007/s00468-023-02448-8}

Lorenz, M.O. (1905) Methods of measuring the concentration of wealth. 
\emph{Journal of the American Statistical Association} 9(70), 209\eqn{-}219. \doi{10.2307/2276207}

Nelder, J.A., Mead, R. (1965) A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Sarabia, J.-M. (1997) A hierarchy of Lorenz curves based on the generalized Tukey's lambda distribution. 
\emph{Econometric Reviews} 16, 305\eqn{-}320. \doi{10.1080/07474939708800389} 

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Sitthiyot, T., Holasut, K. (2023) A universal model for the Lorenz curve with novel applications 
for datasets containing zeros and/or exhibiting extreme inequality. \emph{Scientific Reports} 
13, 4729. \doi{10.1038/s41598-023-31827-x}
}
\seealso{
\code{\link{LeafSizeDist}}, \code{\link{MPerformanceE}}, \code{\link{SarabiaE}}, 
\code{\link{SCSE}}, \code{\link{SHE}}
}
\examples{
  data(LeafSizeDist)

  CulmNumber <- c()
  for(i in 1:length(LeafSizeDist$Code)){
    temp <- as.numeric( strsplit(LeafSizeDist$Code[i], "-", fixed=TRUE)[[1]][1] )
    CulmNumber <- c(CulmNumber, temp)
  }
  uni.CN <- sort( unique(CulmNumber) )  
  ind    <- CulmNumber==uni.CN[1]
  A0     <- LeafSizeDist$A[ind]

  ini.val1 <- list(0.5, 0.1, c(0.01, 0.1, 1, 5, 10), 1, 1)
  ini.val2 <- list(0.5, 0.1, c(0.01, 0.1, 1, 5, 10))
  resu1 <- fitLorenz(MPerformanceE, z=A0, ini.val=ini.val1, simpver=4, fig.opt=TRUE)
  resu2 <- fitLorenz(MPerformanceE, z=A0, ini.val=ini.val2, simpver=5, fig.opt=TRUE)
  resu1$par
  resu2$par

  ini.val3 <- list(0.9, c(10, 50, 100, 500), 1, 0)   
  resu3 <- fitLorenz( SarabiaE, z=A0, ini.val=ini.val3, par.limit=TRUE, 
                      fig.opt=TRUE, control=list(reltol=1e-20, maxit=10000) )
  resu3$par
  resu3$GC

  graphics.off()
}