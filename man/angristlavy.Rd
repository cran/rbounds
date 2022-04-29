\name{angristlavy}
\docType{data}
\alias{angristlavy}
\title{Angrist and Lavy (1999) Maimonides' Rule Data Set}

\description{
  Angrist and Lavy (1999) data set used by Rosenbaum (2010) to
  demonstrate his instrumental variable sensitivity analysis.
}
\usage{data(angristlavy)}
\format{A \code{data.frame} with 172 observations on the following variables).

clasz: Size of class or classes for each cohort.

avgmath: Average math test score for each class.

z: A recode of enrollment with 1 indicating a cohort with 41 or more students, i.e. two classes. 
This serves as the instrument which encourages smaller classes.

pair: matched pair id

}
  
\source{
  Angrist, Joshua and Lavy, Victor (1999). ``Using Maimonides' Rule to
  Estimate the Effect of Class Size on Scholastic
  Achievement.'' \emph{Quarterly Journal of Economics} 114, 533--575.

  Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.
}
\references{
  Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.
}
\keyword{datasets}
