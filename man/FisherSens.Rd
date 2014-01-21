\name{FisherSens}
\alias{FisherSens}

\title{Rosenbaum Sensitivity Analysis for Fisher's Exact Test}
\description{
Calculates sensitivity to hidden bias for Fisher's exact test for a two-by-two contingency table, following the method described in Rosenbaum (2002, sec. 4.4).
}
\usage{
FisherSens(totalN, treatedN, totalSuccesses, treatedSuccesses, Gammas)
}
\arguments{
  \item{totalN}{total number of observations}
  \item{treatedN}{number of treated observations}
  \item{totalSuccesses}{total number of ``successes''}
  \item{treatedSuccesses}{number of successes in treatment group}
  \item{Gammas}{vector of Gammas (bounds on the differential odds of treatment) at which to test the significance of the results}
}

\value{
Returns a matrix with three columns and number of rows equal to the length of "Gammas". Each row indicates the upper and lower bounds for the (one-sided) p-value for a given value of Gamma.
}
\references{Paul R. Rosenbaum. Observational Studies. Springer, New York, 2nd edition, 2002, sec. 4.4}
\author{
Devin Caughey, MIT, \email{caughey@mit.edu}
}

\seealso{See also \code{\link{data.prep}}, \code{\link{binarysens}},
  \code{\link{hlsens}}, \code{\link[Matching]{Match}},
  \code{\link{mcontrol}}
}

\examples{
## Fisher's Lady Tasting Tea: milk first or tea first?
LadyTastingTea <- matrix(c(4, 0, 0, 4), nrow = 2,
                         dimnames = list(Guess = c("Milk", "Tea"),
                                         Truth = c("Milk", "Tea")))
## Define "Milk" as "treated"/"success"
FisherSens(totalN = sum(LadyTastingTea), 
           treatedN = sum(LadyTastingTea["Milk", ]),
           totalSuccesses = sum(LadyTastingTea[, "Milk"]), 
           treatedSuccesses = sum(LadyTastingTea["Milk", "Milk"]), 
           Gammas = seq(1, 2, .2))
## Interpretation: Rejection of the null hypothesis 
## (that the lady cannot discriminate between milk-first and tea-first) 
## is insensitive to bias as large as Gamma = 2.
}

