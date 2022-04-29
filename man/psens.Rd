\name{psens} 
\alias{psens} 
\title{Rosenbaum Bounds for Sign Rank} 

\description{Function to calculate Rosenbaum bounds for continuous or
ordinal outcomes based on Wilcoxon sign rank test.}

\usage{
# Default Method 
psens(x, y, Gamma = 6, GammaInc = 1)
} 

\arguments{ 
  \item{x}{Treatment group outcomes in same order as treatment group.}
  \item{y}{Control group outcomes in same order as treatment group.}
  \item{Gamma}{Upper-bound on gamma parameter.} 
  \item{GammaInc}{To set user-specified increments for gamma parameter.}
}

\references{Rosenbaum, Paul R. (2002) \emph{Observational Studies}.
  Springer-Verlag.} 

\author{Luke Keele, University of Pennsylvania, \email{luke.keele@gmail.com}}
 
\seealso{See also \code{\link{binarysens}},
  \code{\link{hlsens}},
  \code{\link{mcontrol}}
}

\examples{ 
# Replication of Rosenbaum Sensitivity Tests From Chapter 4 of
# Observational Studies

# Data:  Matched Data of Lead Blood Levels in Children
trt <- c(38, 23, 41, 18, 37, 36, 23, 62, 31, 34, 24, 14, 21, 17, 16, 20,
         15, 10, 45, 39, 22, 35, 49, 48, 44, 35, 43, 39, 34, 13, 73, 25,
         27)
ctrl <- c(16, 18, 18, 24, 19, 11, 10, 15, 16, 18, 18, 13, 19, 10, 16,
          16, 24, 13, 9, 14, 21, 19, 7, 18, 19, 12, 11, 22, 25, 16, 13,
          11, 13)

psens(trt, ctrl)
} 
