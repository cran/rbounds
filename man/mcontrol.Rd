\name{mcontrol} 
\alias{mcontrol} 
\title{Sensitivity Analysis For Multiple Matched Controls} 
\description{Function to calculate Rosenbaum bounds for continuous or
ordinal outcomes based on Wilcoxon sign rank test \emph{p}-value when
there are multiple matched control units.}

\usage{
# Default Method 
mcontrol(y, grp.id, treat.id, group.size = 3, Gamma = 4, GammaInc = 1)
} 

\arguments{ 
  \item{y}{Vector of grouped matched outcomes.}
  \item{treat.id}{A vector indicating the treated unit in each matched
    group.}
  \item{grp.id}{A vector indicating matched groups.} 
  \item{group.size}{The size of the matched groups. Three for one
    treated unit and two control units.}
  \item{Gamma}{Upper-bound on gamma parameter.} 
  \item{GammaInc}{To set increments for gamma parameter.} 
} 

\details{The matched data needs to be in a very particular form for
this function to work.  The data must be sorted by matched groups with
indicators for each matched group and for treated and control units.
The simplest way to prepare the data is to use the Match() function and
use the data.prep() function to format the data.

Currently this function only takes matched data with 2 or 3 controls
units matched to each treated unit.

This function does cannot handle data where the number of control units
is not the same for every treated unit.}

\references{Rosenbaum, Paul R. (2002) \emph{Observational Studies}. 
  Springer-Verlag.} 

\author{Luke Keele, University of Pennsylvania, \email{luke.keele@gmail.com}}
 
\seealso{See also \code{\link{binarysens}}, \code{\link{psens}}, 
  \code{\link{hlsens}}}

\examples{ 

grp <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 
         6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10,10,10,11,
         11,11,12,12,12)
trt <- rep(c(1,0,0), 12)
score <- c(-11.39,-8.45, -19.57,-8.33, 3.06, -19.93,-18.73,-11.99,-7.55, 
            11.94, 9.4, -25.16,-0.77, -10.46,-7.27, 24.03, -8.23, 2.67,  
            -4.04, -6.67, -1.12, -14.4, -26.21,5, -1.7, -15.3, -7.73, 
            -0.87, -19.71,-12.69, -3.36, -11.21,-35.83,5.89, -10.79,2)


mcontrol(score, grp, trt, group.size=3)

} 
