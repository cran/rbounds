\name{hlsens} 
\alias{hlsens} 
\title{Rosenbaum Bounds for Hodges-Lehmann Point Estimate} 
\description{ 
Function to calculate Rosenbaum bounds for continuous or ordinal outcomes based on Hodges-Lehmann point estimate. 
} 
\usage{

#Default Method 
hlsens(x, y=NULL, pr=.1, Gamma=6, GammaInc=1)

} 

\arguments{
\item{x}{Treatment group outcomes in same order as treatment group outcomes or an objects from Match().} 
\item{y}{Control group outcomes in same order as treatment group outcomes unnecessary when using Match() object.} 
\item{pr}{Search precision parameter.} 
\item{Gamma}{Upper-bound on gamma parameter.} 
\item{GammaInc}{To set user specified increments for gamma parameter. } 
} 

\details{
For large data sets this function can be quite slow if pr is set to low. If the data set is larger, it is best to set pr to .5 before trying values such as .01.  Generally, the results from the function are insensitive to the value for pr.
}

\references{Rosenbaum, Paul R. (2002) \emph{Observational Studies}. Springer-Verlag.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu}}
 
\seealso{See also \code{\link{data.prep}}, \code{\link{binarysens}},  \code{\link{psens}}, \code{\link[Matching]{Match}}, \code{\link{mcontrol}}}  

\examples{ 
#Replication of Rosenbaum Sensitivity Tests From Chapter 4 of Observational Studies

#Data:  Matched Data of Lead Blood Levels in Children
trt <- c(38,23,41,18,37,36,23,62,31,34,24,14,21,17,16,20,15,10,45,39,22,35,49,48,44,35,43,39,34,13,73,25,27)

ctrl <- c(16,18,18,24,19,11,10,15,16,18,18,13,19,10,16,16,24,13,9,14,21,19,7,18,19,12,11,22,25,16,13,11,13)

hlsens(trt, ctrl)

#
#Load Matching Software and Data
#
library(Matching)
data(lalonde)

#
# Estimate Propensity Score
#
DWglm  <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)

#
#save data objects
#
Y  <- lalonde$re78   #the outcome of interest
Tr <- lalonde$treat #the treatment of interest

#
# Match - without replacement
#             
mDW  <- Match(Y=Y, Tr=Tr, X=DWglm$fitted, replace=FALSE)

#
# One should check balance, but let's skip that step for now.
#

#
# Sensitivity Test
#
hlsens(mDW, pr=.1, Gamma=2, GammaInc=.1)
} 
