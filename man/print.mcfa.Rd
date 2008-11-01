\name{print.mcfa}
\alias{print.mcfa}
\title{Print an object of the class mcfa}
\description{
  Printing method for one of two possible objects returned by  \code{cfa()}
}
\usage{
\method{print}{mcfa}(x,...)
}
\arguments{
  \item{x}{An object of the class mcfa}
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

\note{Note that \code{cfa()} will return an object with the class scfa if there is only one
row of counts. If there are two or more of them, an object with the class mcfa is returned.
In contrast \code{scfa()} and \code{mcfa()} return a list which has no class of it's own.}

\seealso{ \code{\link{cfa}}, \code{\link{mcfa}}}

\examples{
# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts1<-trunc(runif(250)*10) 
counts2<-trunc(runif(250)*10)
result<-cfa(configs,cbind(counts1,counts2))
print(result)
}

\keyword{multivariate}
\keyword{htest}


