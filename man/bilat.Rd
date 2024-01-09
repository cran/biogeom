\name{bilat}
\encoding{UTF-8}
\alias{bilat}
\title{
Measure of the Extent of Bilateral Symmetry of A Polygon
}
\description{
\code{bilat} is used to measure the extent of bilateral (a)symmetry and other measures for a polygon (e.g., a leaf).
}
\usage{
bilat(x, y, strip.num = 200, peri.np = NULL, n.loop = 60,
      auto.search = TRUE, animation.fig = TRUE, time.interval = 0.001,   
      unit = "cm", main = NULL, diff.fig = TRUE, angle = NULL,  
      ratiox = 0.02, ratioy = 0.08, fd.opt = TRUE, frac.fig = TRUE,
      denomi.range = seq(8, 30, by = 1))
}
\arguments{
  \item{x}{the \eqn{x} coordinates of a polygon's boundary.}
  \item{y}{the \eqn{y} coordinates of a polygon's boundary.}
  \item{strip.num}{the number of equidistant strips intersecting with 
the polygon that are horizontally placed. See Shi et al. (2018, 2020) for details.}  
  \item{peri.np}{the number of data points on the boundary retained for calculating the perimeter of the polygon.}
  \item{n.loop}{the number of data points to randomly sample for calculating the mean perimeter of the polygon.}
  \item{auto.search}{an optional argument to automatically search the maximum distance between two points on the polygon's boundary.}
  \item{animation.fig}{the option of showing the data points on the polygon's boundary in an animation.}
  \item{time.interval}{the time interval at which to suspend execution, in seconds.}
  \item{unit}{the units of the \eqn{x}-axis and the \eqn{y}-axis when showing the polygon.}
  \item{main}{the main title of the figure.}  
  \item{diff.fig}{an optional argument to draw the differences in areas between the intersections of the strips 
with the upper part of the polygon and the intersections of the strips with the lower part of the polygon. 
The polygon is divided into the upper and lower parts by the \eqn{x}-axis. See Shi et al. (2018, 2020) for details.}
  \item{angle}{the angle between the major axis (i.e., the leaf length axis) and the \eqn{x}-axis, which can be defined by the user.}
  \item{fd.opt}{An optional argument to use the box-counting method to calculate the fractal
dimension of the polygon's boundary on a log-log scale.}
  \item{ratiox}{the \eqn{x} coordinate of the location parameter for positioning the legend in the plot of the linear fitting.}
  \item{ratioy}{the \eqn{y} coordinate of the location parameter for positioning the legend in the plot of the linear fitting.}
  \item{frac.fig}{an optional argument to draw the results of the linear fitting using the box-counting method to calculate the fractal
dimension of the polygon's boundary on a log-log scale.}
  \item{denomi.range}{the number of equidistant segments of the maximum range 
between the range of the \eqn{x} coordinates and that of the \eqn{y} coordinates.}

}
\details{
The data of \code{x} and \code{y} should be the coordinates adjusted using the \code{\link{adjdata}} function. 
If \code{peri.np = NULL}, the number of \code{length(x)} is used to calculate the perimeter of the polygon; 
if \code{peri.np} is a positive integer, the number of data points retained on the polygon's boundary  
is equal to \code{peri.np} and random sampling for retaining \code{peri.np} data points is carried out 
\code{n.loop} times for calculating the mean perimeter of the polygon. That is to say, the final output for 
the perimeter is the mean of the \code{n.loop} perimeters (i.e., replicates). If the user wants to get a consistent result 
for the mean perimeter, the \code{\link{set.seed}} function can be used. In addition, if \code{length(x) < peri.np}, 
\code{peri.np} then becomes \code{length(x)} rather than the specified value in Arguments. 
If the polygon apparently has a major axis (e.g., the leaf length axis for an ovate leaf), \code{auto.search} 
is appropriate. If the major axis of the polygon is not the straight line through two points on the polygon's 
boundary having the maximum distance, the user can define the major axis using the \code{\link{locator}} 
function in \pkg{graphic} by clicking two points on or near the polygon's boundary. The location of the first click 
should be northeast of the location of the second click. This means that the angle between the straight line 
through the locations of the two clicks and the \eqn{x}-axis should range from 0 to \eqn{\pi}/2. The locations of the 
clicks can be on the boundary or be approximate to the boundary. The function will automatically find the nearest 
data point on the boundary to the location of each click. 
When \code{angle = NULL}, the observed polygon will be shown at its initial angle in the scanned image; 
when \code{angle} is a numerical value (e.g., \eqn{\pi/4}) defined by the user, it indicates that the major axis 
is rotated \eqn{\pi/4} counterclockwise from the \eqn{x}-axis. 
}
\value{
  \item{x}{the \eqn{x} coordinates retained on the polygon's boundary.}
  \item{y}{the \eqn{y} coordinates retained on the polygon's boundary.}
   \item{phi}{the angle between the length axis (i.e., the major axis) of the polygon and the \eqn{x}-axis.}
  \item{n1}{the number of data points on the upper boundary of the polygon.}
  \item{n2}{the number of data points on the lower boundary of the polygon.}
  \item{n}{the number of data points on the whole polygon's boundary.}
  \item{total.poly}{an object of class "ppp" representing a point pattern dataset 
in the two-dimensional plane, representing the polygon's boundary.}
  \item{upper.poly}{an object of class "ppp" representing a point pattern dataset 
in the two-dimensional plane, representing the upper boundary of the polygon along the \eqn{x}-axis.}
  \item{lower.poly}{an object of class "ppp" representing a point pattern dataset 
in the two-dimensional plane, representing the lower boundary of the polygon along the \eqn{x}-axis.}
  \item{D}{the differences in areas between the upper and lower boundaries of the polygon.}
  \item{par.upper.area}{the area of the upper boundary of the polygon along the \eqn{x}-axis.}
  \item{par.lower.area}{the area of the lower boundary of the polygon along the \eqn{x}-axis.}
  \item{SI}{the standardized index for bilateral (a)symmetry for the polygon.}
  \item{AR}{the ratio of the areas of the upper to the lower parts of the polygon.}
  \item{scan.length}{the length of the polygon. The default is the maximum distance between two points on 
the polygon's boundary.}
  \item{scan.width}{the maximum width of the polygon.}
  \item{scan.area}{the area of the polygon.}
  \item{scan.perimeter}{the perimeter of the polygon 
based on all data points or a mean of \code{n.loop} 
replicates of perimeters using the \code{peri.np} data points retained on the polygon's boundary.}
  \item{x.width}{distance from the base to a point on the major axis 
associated with the maximum width of the polygon.}
  \item{width.1e}{the width associated with 1/8 of \code{scan.length} (starting from the base of the polygon).} 
  \item{width.2e}{the width associated with 2/8 of \code{scan.length} (starting from the base of the polygon).}
  \item{width.4e}{the width associated with 4/8 of \code{scan.length} (starting from the base of the polygon).}
  \item{width.6e}{the width associated with 6/8 of \code{scan.length} (starting from the base of the polygon).}
  \item{width.7e}{the width associated with 7/8 of \code{scan.length} (starting from the base of the polygon).}
  \item{bi.test}{the testing results for \code{D} using the Wilcoxon signed rank test with continuity correction. 
See the \code{\link{wilcox.test}} function in \pkg{stats}.}
  \item{a}{the estimate of the intercept obtained using the box-counting method to calculate the fractal dimension of the polygon's boundary.}
  \item{sd.a}{the standard deviation of the estimated intercept.}
  \item{lci.a}{the lower bound of the 95\% confidence interval of the estimated intercept.}
  \item{uci.a}{the upper bound of the 95\% confidence interval of the estimated intercept.}
  \item{b}{the estimate of the slope obtained using the box-counting method to calculate the fractal dimension of the polygon's boundary.}
  \item{sd.b}{the standard deviation of the estimated slope.}
  \item{lci.a}{the lower bound of the 95\% confidence interval of the estimated slope.}
  \item{uci.a}{the upper bound of the 95\% confidence interval of the estimated slope.}
  \item{r.sq}{the coefficient of determination obtained when using the box-counting method to calculate the fractal dimension of the polygon's boundary.}
  \item{delta}{the vector of box sizes used in the box-counting method to calculate the fractal dimension of the polygon's boundary.}
  \item{N}{the number of boxes that include at least one pixel of the polygon's boundary.}
}
\note{
The polygon is expected to have an apparent major axis (e.g., the straight line through two points 
on the polygon's boundary having the maximum distance or one that can be clearly defined to pass by 
two landmarks on the polygon's boundary [i.e., the leaf length axis, the egg length axis, etc.]). 
The polygon is placed with its major axis overlapping the \eqn{x}-axis; the base of the polygon is 
located at the origin; the apex of the polygon is located to the right of the base. 
\code{phi} is equal to \code{angle} when \code{angle} is not null. 
In theory, \code{n1 + n2 = n}, but in most cases \code{n1 + n2} is slightly smaller than \code{n}. 
The reason is that very few boundary points fall outside the the lower and upper boundaries of the polygon 
when using the \code{\link{intersect.owin}} function in \pkg{spatstat.geom}. However, this does not considerably affect the results. 
The log-transformed \code{SI} and 
the log-transformed \code{AR} are demontrated to have a more symmetrical frequency distribution than their original forms. 
This is important when performing an analysis of variance between (or among) groups to compared 
their extents of bilateral (a)symmetry. See Shi et al. (2020) for details. 
The box-counting approach uses a group of boxes (squares for simplicity) with different
sizes (\eqn{\delta}) to divide the leaf vein image into different parts. Let \eqn{N} represent the number
of boxes that include at least one pixel of the polygon's boundary. 
The maximum of the range of the \eqn{x} coordinates and the range of the \eqn{y} coordinates 
for the pixels of the polygon's boundary is defined as \eqn{z}. Let \eqn{\delta} represent the vector of 
\eqn{z}/\code{denomi.range}. We then used the following equation to calculate the fractal 
dimension of the polygon's boundary: 
\deqn{\mathrm{ln } N = a + b\, \mathrm{ln} \left(\delta^{-1}\right),}
where \eqn{b} is the theoretical value of the fractal dimension. We can use its estimate as the 
numerical value of the fractal dimension for the polygon's boundary. 
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Gielis, J., Quinn, B.K., Niklas, K.J., Ratkowsky, D.A., Schrader, J., Ruan, H., 
Wang, L., Niinemets, Ü. (2022) 'biogeom': An R package for simulating and fitting natural 
shapes. \emph{Annals of the New York Academy of Sciences} 1516, 123\eqn{-}134. \doi{10.1111/nyas.14862}

Shi, P., Niinemets, Ü., Hui, C., Niklas, K.J., Yu, X., Hölscher, D. (2020) 
Leaf bilateral symmetry and the scaling of the perimeter vs. the surface area in 15 vine species. 
\emph{Forests} 11, 246. \doi{10.3390/f11020246}

Shi, P., Zheng, X., Ratkowsky, D.A., Li, Y., Wang, P., Cheng, L. (2018)  
A simple method for measuring the bilateral symmetry of leaves. 
\emph{Symmetry} 10, 118. \doi{10.3390/sym10040118}
}
\seealso{
\code{\link{adjdata}}, \code{\link{fracdim}}
}
\examples{
data(bambooleaves)

uni.C <- sort( unique(bambooleaves$Code) )
ind   <- 3
Data  <- bambooleaves[bambooleaves$Code==uni.C[ind], ]
x0    <- Data$x
y0    <- Data$y

dev.new()
plot( x0, y0, asp=1, type="l", cex.lab=1.5, cex.axis=1.5, 
      xlab=expression(italic(x)), ylab=expression(italic(y)) )

Res1 <- adjdata(x0, y0, ub.np=2000, len.pro=1/20)
x1   <- Res1$x
y1   <- Res1$y


\donttest{
  Res2 <- bilat( x=x1, y=y1, time.interval=0.00045, 
                 peri.np=NULL, auto.search=TRUE, 
                 fd.opt=TRUE )
  Res2$scan.perimeter

  set.seed(123)
  Res3 <- bilat( x=x1, y=y1, time.interval=0.00045, 
                 peri.np=500, n.loop=30,
                 auto.search=TRUE, fd.opt=FALSE )
  Res3$scan.perimeter

  set.seed(123)
  Res4 <- bilat( x=x1, y=y1, time.interval=0.00045, 
                 peri.np=500, n.loop=30,
                 auto.search=TRUE, fd.opt=FALSE, angle=pi/4 )
  Res4$scan.perimeter

  set.seed(123)  
  Res5 <- bilat( x=x1, y=y1, time.interval=0.00045, 
                 peri.np=500, n.loop=30,
                 auto.search=TRUE, fd.opt=FALSE, angle=0 )
  Res5$scan.perimeter
  
  if(interactive()){
    # The angle between the leaf length axis (namely the straight 
    # line through the leaf apex and base) and the horizontal axis 
    # should be between 0 and pi/2 for a scanned leaf's profile. 
    # Here, the user needs to first click the leaf apex, 
    # and then click the leaf base. 
    set.seed(123)  
    Res6 <- bilat( x=x1, y=y1, time.interval=0.00045, 
                 peri.np=500, n.loop=30,
                 auto.search=FALSE, fd.opt=FALSE, angle=NULL )
    Res6$scan.perimeter
  }

  set.seed(NULL)
}

graphics.off()
}