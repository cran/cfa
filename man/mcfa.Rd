\name{mcfa}
\alias{mcfa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Repeated measures analysis of configuration frequencies }
\description{
Determines the frequency of all combinations of variable values (their 
configurations) over time or with and without treatment in comparison with their expected frequency and displays them in the order of decreasing chi-squared. In addition, a global chi squared is calculated.
}
\usage{
mcfa(configmatrix, cntmatrix, descending=TRUE, sortonchisq=TRUE, ignore.na=FALSE, verbose=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{configmatrix}{Dataframe with the variables to be analyzed}
 \item{cntmatrix}{Matrix with >=2 columns of counts (containing 1 if the data are not aggregated) }
 \item{descending}{Output in the order of decreasing chi squared }
 \item{sortonchisq}{Sort output on chi squared}  
 \item{ignore.na}{Ignore (casewise) missing data in the configurations}
 \item{verbose}{Long output}
}
\details{
Each variable must have at least two different values and may have more 
(extension of the classical CFA). The \code{configmatrix} must consist of at least 
two variables (columns). Factors and numbers are both accepted (the numbers 
are internally converted to factors). The \code{cntmatrix} must be 
numeric.
Counts should be at least = 5 for the chi squared test to be reliable but when 
using the CFA as a purely heuristic tool counts of 0 are possible. 
}
\value{
 A list with class \code{"mcfa"} contains the tabular results and the overall parameters
  \item{Row names}{Configuration}
  \item{n(1)..n(configs)}{Frequency (count) of this configuration}
  \item{expected(1)..expected(config)}{Expected Frequency (count) of this configuration}
  \item{chi.sq}{Chi squared for the given configuration}
  \item{p(chisq)}{p(chi squared) for the given configuration}
  \item{Overall chi squared }{Overall chi squared }
  \item{p(chi squared)  }{ p(overall chi squared) }
  \item{Degrees of freedom }{ Overall degrees of freedom }
}

\references{ 
Krauth J., Lienert G. A. (1973, Reprint 1995)
Die Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in Psychologie und Medizin
Beltz Psychologie Verlagsunion

Eye, A. von (1990)
Introduction to configural frequency analysis. The search for types and anti-types in cross-classification. 
Cambride 1990
}
\author{Stefan Funke <funke@attglobal.net>}


\section{WARNING }{The program is implemented in R itself rather than a 
compiled library and therefore slow. In most cases the input is a 
pre-aggregated table and speed is no problem because the \code{configmatrix} 
is small. There are no hard-coded limits in the program so even large tables 
can be processed but this will take time and memory. 

The outout table can be very wide if the levels of the factors variables are 
long so \code{options(width=..)} may need to be adjusted}

\seealso{ \code{\link{cfa}}, \code{\link{hier.cfa}}, \code{\link{boot.cfa}}
}
 
\examples{
library(cfa)
data(cfa2dat)
mcfa(cfa2dat[1:3],cfa2dat[4:5],verbose=T)
}
\keyword{ mcfa}%-- one or more ...
