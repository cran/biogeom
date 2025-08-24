\name{fitGS}
\encoding{UTF-8}
\alias{fitGS}
\title{
Data-Fitting Function for the Geometric Series of Size Distribution Measurements
}
\description{
\code{fitGS} is used to estimate the model parameters of a geometric series, i.e., the first term and common ratio.
}
\usage{
fitGS(A, ini.val = NULL, control = list(), par.list = FALSE, fig.opt = TRUE)
}
\arguments{
  \item{A}{A sequence of size distribution measurements, e.g., the temporal or spatial progression of leaf area in an individual plant.}
  \item{ini.val}{the list of initial values for the model parameters.}
  \item{control}{the list of control parameters for using the \code{\link[stats]{optim}} function in package \pkg{stats}.}
  \item{par.list}{an optional argument to show the list of parameters on the screen.}
  \item{fig.opt}{an optional argument to draw (i) the observed and predicted \eqn{y} values 
(size distribution measurements) sorted in ascending order, and (ii) the the observed and predicted \eqn{z} values 
(cumulative size distribution measurements) sorted in ascending order.}
}
\details{
In general, there is no need to set the initial values for model parameters by users (i.e., \code{ini.val = NULL}). 
The approach proposed by Yan et al. (2025) is used to find the suitable initial values for model parameters. 
The Nelder-Mead algorithm (Nelder and Mead, 1965) is used to carry out the optimization of minimizing the 
residual sum of squares (RSS) between the observed and predicted \eqn{y} values (Deng et al., 2025; Yan et al. 2025). 
The \code{\link[stats]{optim}} function in package \pkg{stats} was used to carry out the Nelder-Mead algorithm. 
\code{fig.opt = TRUE} generates two panels (i.e., \eqn{y} vs. \eqn{\hat{y}}, and \eqn{z} vs. \eqn{\hat{z}}) in a figure.
}
\value{
  \item{x}{the ascending number from 1 to \code{length(A)}.}
  \item{y}{the observed size distribution measurements (\eqn{y} values) sorted in ascending order.}
  \item{y.theo}{the predicted size distribution measurements (\eqn{y} values) sorted in ascending order by the geometric seires.}
  \item{z}{the observed cumulative size distribution measurements (\eqn{z} values) sorted in ascending order.}
  \item{z.theo}{the predicted cumulative size distribution measurements (\eqn{z} values) sorted in ascending order by the geometric seires.}
  \item{par}{the estimates of the model parameters corresponding to the first term and common ratio, respectively.}
  \item{r.sq}{the coefficient of determination between the observed and predicted \eqn{y} values.}
  \item{RSS}{the residual sum of squares between the observed and predicted \eqn{y} values.}
  \item{sample.size}{the number of data points used in the data fitting.}
  \item{RMSE1}{the root mean square error between the observed and predicted \eqn{y} values.}
  \item{MAPE1}{the mean absolute percent error between the observed and predicted \eqn{y} values (in \%).}
  \item{RMSE2}{the root mean square error between the observed and predicted \eqn{z} values.}
  \item{MAPE2}{the mean absolute percent error between the observed and predicted \eqn{z} values (in \%).}
}
\note{
In the outputs, the parameter vector \code{par} was estimated based on the observed 
and predicted \eqn{y} values rather than the observed and predicted \eqn{z} values.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{

Deng, L., Wang, J., Zhang, L., Hölscher, D., Shi, P. (2025) Testing the validity of 
the Koyama-Smith equation and the power-law equation using 3231 tepals of a \emph{Magnolia} 
species. \emph{Trees \eqn{-} Structure and Function} 39, 74. \doi{10.1007/s00468-025-02645-7}

Nelder, J.A., Mead, R. (1965) A simplex method for function minimization. 
\emph{Computer Journal} 7, 308\eqn{-}313. \doi{10.1093/comjnl/7.4.308}

Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Yan, C., Shi, P., Yao, W., Yu, K., Niinemets, Ü. (2025) A nonlinear fitting method provides 
strong support for geometric series of stomatal area in 12 Magnoliaceae species. 
\emph{Plants} 14, 893. \doi{10.3390/plants14060893}
}
\examples{
  # A sequence of tepal area measurements of a Magnolia flower (Deng et al., 2025)
  A   <- c(56.65, 43.37, 49.61, 56.27, 56.66, 49.45, 46.56, 43.42, 44.80)
  ReS <- fitGS(A)
  names( ReS)

  graphics.off()
}