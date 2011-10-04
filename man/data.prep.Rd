\name{data.prep} 
\alias{data.prep} 
\title{Format Match() object for sensitivity test.} 
\description{This function reshapes the output from Match() to create
  the necessary objects for mcontrol().}
\usage{
# Default Method 
data.prep(obj, Y = NULL, group.size = 3)
} 

\arguments{ 
  \item{obj}{An object from the Match() function.}
  \item{Y}{Internal argument to the Match object.}
  \item{group.size}{The size of the matched groups. Three for one
    treated unit and two control units.}
}

\details{This functions takes a Match() object and formats it for use
with the mcontrol() function.  The output is a list with the three
objects needed for the arguments of the mcontrol() function.}

\value{
  \item{Y}{The matched outcomes}
  \item{id}{A vector which identifies the matched groups: 1, 1, 1 for
    matched group one; 2, 2, 2 for match group 2, etc.}  
  \item{treat}{A vector with 1's for treated units and 0's for control
    units}
}


\references{Rosenbaum, Paul R. (2002) \emph{Observational Studies}.
  Springer-Verlag.} 

\author{Luke Keele, Penn State University, \email{ljk20@psu.edu}} 

\seealso{See also \code{\link{binarysens}}, \code{\link{psens}},
  \code{\link{hlsens}}, \code{\link[Matching]{Match}},
  \code{\link{mcontrol}}} 

\examples{ 

# Load Matching Software and Data
library(Matching)
data(lalonde)

# Estimate Propensity Score
DWglm <- glm(treat~age + I(age^2) + educ + I(educ^2) + black + hisp +
             married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)

# Save data objects
Y  <- lalonde$re78    # the outcome of interest
Tr <- lalonde$treat   # the treatment of interest

# Match      
mDW  <- Match(Y = Y, Tr = Tr, X = DWglm$fitted, M = 2)

# One should check balance, but let's skip that step for now.

# Create Data Object
tmp <- data.prep(mDW, group.size=3)

# Sensitivity Test
mcontrol(tmp$Y, tmp$id, tmp$treat, group.size = 3)
} 
