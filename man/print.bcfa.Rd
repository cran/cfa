\name{print.bcfa}
\alias{print.bcfa}
\title{Print an object of the class hcfa}
\description{
  Printing method for an object returned by  \code{boot.cfa()}
}
\usage{
\method{print}{bcfa}(x,...)
}
\arguments{
  \item{x}{An object of the class bcfa}
  \item{...}{Additional arguments given to \code{print}}
}
\details{
  This function is usually called implicitely.
}
\value{Returns NULL
}

\references{Krauth J., Lienert G. A. (1973, Reprint 1995) Die
     Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in
     Psychologie und Medizin, Beltz Psychologie Verlagsunion
     
     Lautsch, E., von Weber S. (1995) Methoden und Anwendungen der Konfigurationsfrequenzanalyse 
     in Psychologie und Medizin, Beltz Psychologie Verlagsunion

     Eye, A. von (1990) Introduction to configural frequency analysis.
     The search for types and anti-types in cross-classification. 
     Cambride 1990 }


\author{Stefan Funke <s.funke@t-online.de>}

\seealso{ \code{\link{bcfa}}}

\examples{
# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts<-trunc(runif(250)*10)
result<-bcfa(configs,counts,runs=25) 
print(result)
}

\keyword{multivariate}
\keyword{htest}


