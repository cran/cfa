\name{hier.cfa}
\alias{hier.cfa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Hierachical analysis of configuration frequencies}
\description{Recursively eliminates variables to generate subtables which are subjected 
to a CFA and calculates the overall chi squared.
}
\usage{
hier.cfa(configmatrix, cntvector)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{configmatrix}{ Dataframe with the variables to be analyzed}
 \item{cntvector}{ Vector of counts (1 if the data are not aggregated)}
}
\details{The hierarchical CFA is used to eliminate unnecessary variables. 
The higher the overall chi squared, the better the information from the 
original table is preserved. Given a constant number of variables in the 
subtables, the subtable with the highest chi squared is the "best". 
Sometimes it makes sense to accept a somewhat lower chi squared if the number 
of variables is lower because of better interpretability.
}
\value{
  \item{Row names}{Configuration}
  \item{Overall chi squared}{Overall chi squared for this subtable}
  \item{order }{Number of variables used}

}
\references{
Krauth J., Lienert G. A. (1973, Reprint 1995)
Die Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in Psychologie und Medizin
Beltz Psychologie Verlagsunion

Eye, A. von (1990)
Introduction to configural frequency analysis. The search for types and anti-types in cross-classification. 
Cambride 1990
}
\author{Stefan Funke <funke@attglobal.net> }

\section{WARNING }{The program is implemented in R itself rather than a 
library and therefore slow. It repeatedly calls a simplified (and faster) 
version of \code{cfa}. With increasing number of variables in the 
configuration the running time will incease expontentially} 

\seealso{ \code{\link{cfa}}, \code{\link{boot.cfa}} }



\examples{
library(cfa)
data(cfadat)
hier.cfa(cfadat[c("gender","married","children")],cfadat["count"])
}
\keyword{hier.cfa}%-- one or more ...
