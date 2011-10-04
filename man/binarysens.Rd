


\name{binarysens} 
\alias{binarysens} 
\title{Binary Sensitivity Test} 
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

\author{Luke Keele, Penn State University, \email{ljk20@psu.edu}}
 
\seealso{See also \code{\link{data.prep}}, \code{\link{psens}},  \code{\link{hlsens}}, \code{\link[Matching]{Match}}, \code{\link{mcontrol}}
}

\examples{ 
# Example From Rosenbaum Observational Studies Pg 112
# Success: Died From Lung Cancer
# 110 Discrepant Pairs
# 12 Discrepant Pairs

# Sensitivity Test
binarysens(12,110)

# Example Using Match()

#Load Matching Software and Data

library(Matching)
data(GerberGreenImai)

# Estimate Propensity Score
pscore.glm <- glm(PHN.C1 ~ PERSONS + VOTE96.1 + NEW + 
    MAJORPTY + AGE + WARD + PERSONS:VOTE96.1 + PERSONS:NEW + 
    AGE2, family = binomial(logit), data = GerberGreenImai)

#save data objects
D <- GerberGreenImai$PHN.C1
Y <- GerberGreenImai$VOTED98
X <- fitted(pscore.glm)

# Match - without replacement
m.obj <- Match(Y = Y, Tr = D, X = X, M = 1, replace=FALSE)
summary(m.obj)

# One should check balance, but let's skip that step for now.

# Sensitivity Test
binarysens(m.obj, Gamma=2, GammaInc=.1)

} 
