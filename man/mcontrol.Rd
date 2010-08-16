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

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} }
 
\seealso{See also \code{\link{data.prep}}, \code{\link{binarysens}},
  \code{\link{psens}}, \code{\link[Matching]{Match}},
  \code{\link{hlsens}}}

\examples{ 

# Load Matching Software and Data
library(Matching)
data(lalonde)

# Estimate Propensity Score
DWglm <- glm(treat~age + I(age^2) + educ + I(educ^2) + black + hisp +
             married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)

# Save data objects
Y  <- lalonde$re78   # the outcome of interest
Tr <- lalonde$treat  # the treatment of interest

# Match - without replacement
mDW  <- Match(Y=Y, Tr=Tr, X=DWglm$fitted, M=2, replace=FALSE)

# One should check balance, but let's skip that step for now.

# Create Data Object:
tmp <- data.prep(mDW, group.size=3)

# Sensitivity Test
mcontrol(tmp$Y, tmp$id, tmp$treat, group.size=3)
} 
