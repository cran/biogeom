\name{PlanCoor}
\encoding{UTF-8}
\alias{PlanCoor}
\title{
Extracts the Planar Coordinates the 2-D Profiles of Apical Meristems
}
\description{
\code{PlanCoor} is used to extract the planar coordinates the 2-D profiles 
of shoot or root apical meristems in grayscale bitmap (.bmp) images saved in a folder.
}
\usage{
PlanCoor( folder.name, lower.val = 0, upper.val = 250, ratio = 1,  
          fig.opt = TRUE, np = NULL, unit = "unitless", verbose = TRUE )
}
\arguments{
  \item{folder.name}{the name of a folder that saves the 2-D profiles 
of shoot or root apical meristems in grayscale .bmp images}
  \item{lower.val}{the lower bound of pixel values for selecting an AM's profile curve from a grayscale image.}
  \item{upper.val}{the upper bound of pixel values for selecting an AM's profile curve from a grayscale image.}
  \item{ratio}{the scaling factor between the actual image dimensions and the desired coordinate system scale.}
  \item{fig.opt}{an optional argument to draw the extracted the AM's profile curve.}
  \item{np}{the number of data points forming an AM's profile curve.}
  \item{unit}{the unit of the \eqn{x}-axis and the \eqn{y}-axis when showing an AM's profile curve.}
  \item{verbose}{an optional argument allowing the user to suppress the printing of computation progress.}
}
\details{
The apical meristem profiles are recommended to be extracted using the Pen Tool in Adobe Photoshop 
(version >= 22.4.2; Adobe, San Jose, CA, USA) and saved as grayscale bitmap (.bmp) images. 
In each image, the profile should be depicted in black or a gray value close to black, 
whereas the remaining areas should be set to white (with pixel values of 255). 
All profiles are recommended to be adjusted to form downward-opening curves using Adobe Photoshop. 
The resulting profiles are to be saved in a folder within the current R working directory. 
The name of this folder is used as an argument \code{folder.name} in the \code{\link{PlanCoor}} function. 
Users can adjust the values of \code{lower.val} and \code{upper.val} to enhance 
the clarity of the shoot or root apical meristem profiles. The argument \code{ratio} 
represents the scaling factor between the actual image dimensions and the desired coordinate system scale. 
For example, if an image has a width of 7 cm and a height of 10 cm, setting \code{ratio = 2} will result in 
\eqn{x}- and \eqn{y}-coordinates being mapped to a window of 14 cm by 20 cm. The default of \code{np} is 
\code{NULL}, and the number of data points forming an AM's profile curve depends on the resolution of the 
image and the size of the extracted AM's profile. If the use defines an integer such as 500, the final 
number of data points is equal to the minimum value of 500 and the number of data point when \code{np = NULL}.
}
\value{
  \item{FileName}{the names of the files in the folder but removing the file extension 
(e.g., "myimage" from a file named "myimage.bmp" by dropping its extension ".bmp").}
  \item{x}{the \eqn{x}-coordinates of the profile curve(s).}
  \item{y}{the \eqn{y}-coordinates of the profile curve(s).}   
}
\note{
In the folder, one image is also allowed. Assume that all grayscale .bmp images are saved in 
a folder named "Image" in the current working directory. Using \code{PlanCoor(folder.name = "Image", 
fig.opt = FALSE)} can obtain the \eqn{x}- and \eqn{y}-coordinates of the profiles. 

\eqn{\quad} The \code{\link{PlanCoor}} function can also be used to extract two-dimensional profiles 
of other objects, such as leaves. The leaf profile can first be obtained using image processing 
software (e.g., Adobe Photoshop). It is recommended to fill the exterior of the leaf profile 
in black and the interior in white. In such cases, the values of \code{lower.val} and \code{upper.val} 
should be adjusted accordingly to ensure a clear leaf profile is extracted.
}
\author{
Peijian Shi \email{pjshi@njfu.edu.cn}, Johan Gielis \email{johan.gielis@uantwerpen.be}, 
Brady K. Quinn \email{Brady.Quinn@dfo-mpo.gc.ca}.
}
\references{
Shi, P., Chen, L., Quinn, B.K., Yu, K., Miao, Q., Guo, X., Lian, M., Gielis, J., Niklas, K.J. (2023) 
A simple way to calculate the volume and surface area of avian eggs. 
\emph{Annals of the New York Academy of Sciences} 1524, 118\eqn{-}131. \doi{10.1111/nyas.15000}
}
\seealso{
\code{\link{fitAM}}, \code{\link{SAMs}}, \code{\link{SurfaceAreaAM}}, \code{\link{VolumeAM}}
}
\examples{
  folder_path1 <- system.file("extdata", "SAM", package="biogeom")
  bmp.res1     <- PlanCoor( folder_path1, lower.val = 0, upper.val = 250, 
                            fig.opt = TRUE, np = NULL )

  folder_path2 <- system.file("extdata", "Leaf", package="biogeom")
  bmp.res2     <- PlanCoor( folder_path2, lower.val = 1, upper.val = 254, 
                            fig.opt = TRUE, np = NULL, unit = "cm", verbose = FALSE )

  graphics.off()
}
