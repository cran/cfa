\name{boot.cfa}
\alias{boot.cfa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Bootstrap CFA}
\description{
Check if the configurations detected are stable. This is done by counting the 
configurations found to be significant in the CFAs performed.
}
\usage{
boot.cfa(configmatrix, cntvector, runs=100, bonferroni=F, sig.limit=0.05)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{configmatrix}{ Dataframe with the variables to be analyzed}
 \item{cntvector}{ Vector of counts (1 if the data are not aggregated}
 \item{runs}{ Number of runs and therefore samples to be generated.}
 \item{bonferroni}{ Use Bonferroni-adjustment in determining if the configuration is considered significant}
 \item{sig.limit}{ Limit of p which cause the configuration to be counted as "signifcant"}
}
\details{
  The repeated use of tests of significance is obviously problematic. The result is therfore to be interpreted in a strictly heuristic way.
}
\value{
 \item{Row names}{Configuration}
 \item{cnt.antitype}{Number of cases where n<expected}
 \item{cnt.type}{Number of cases where n>expected}
 \item{pct.types}{Percentage of cases where n>expected}
 \item{cnt.sig}{Number of cases where the p(z) for the configuration was considered significant}
 \item{pct.cnt.sig}{Percentage of cases where the p(z) for the configuration was considered significant}
 

}
\references{
Lautsch E., von Weber S. (1995)
Methoden und Anwendungen der Konfigurationsfrequenzanalyse (KFA)
Beltz Verlagsunion (1995)
}
\author{Stefan Funke <funke@attglobal.net>}

\section{WARNING }{The program is implemented in R itself rather than a 
compiled library and therefore slow. It repeatedly calls \code{cfa}. 
In addition the number of runs must be high to give a reliable result. 
The default of 100 is the absolute mininum. So a run of the program under 
realistic conditions may take from minutes to hours.
} 

\seealso{ \code{\link{cfa}}, \code{\link{hier.cfa}} }



\examples{
library(cfa)
data(cfadat)
boot.cfa(cfadat[c("gender","married","children")],cfadat["count"])
}
\keyword{ boot.cfa}%-- one or more ...
