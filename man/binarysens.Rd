


\name{binarysens} 
\alias{binarysens} 
\title{Sensitivity Test for Matched Binary Outcome} 
\description{ 
Function to calculate Rosenbaum bounds for binary data. 
} 
\usage{ 
binarysens(x,y, Gamma=6, GammaInc=1)
} 
\arguments{ 
\item{x}{Count of the first set of discrepant pairs in a table of treated and control outcomes.} 
\item{y}{Count of the second set of discrepant pairs in a table of treated and control outcomes.} 
\item{Gamma}{Upper-bound on gamma parameter.} 
\item{GammaInc}{To set user specified increments for gamma parameter. } 
} 

\references{Rosenbaum, Paul R. (2002) \emph{Observational Studies}. Springer-Verlag.} 

\author{Luke Keele, University of Pennsylvania, \email{luke.keele@gmail.com}}
 
\seealso{See also \code{\link{psens}},  \code{\link{hlsens}}, \code{\link{mcontrol}}
}

\examples{ 
# Example From Rosenbaum Observational Studies Pg 112
# Success: Died From Lung Cancer
# 110 Discrepant Pairs
# 12 Discrepant Pairs

# Sensitivity Test
binarysens(12,110)

} 
