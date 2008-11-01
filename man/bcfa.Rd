\name{bcfa}
\alias{bcfa}
\title{Bootstrap-CFA}
\description{ The bootstrap-CFA tries to replicate the pattern of significant configurations by
re-sampling. 
}
\usage{
bcfa(configs, cnts, runs=100, sig.item="sig.z",...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{configs}{ Contains the configurations. This can be a dataframe or a matrix. The dataframe can contain numbers,
              characters, factors, or booleans. The matrix can consist of numbers, characters or booleans 
              (factors are implicitely re-converted to numerical levels). There must be >=3 columns.
            }
  \item{cnts}{ Contains the counts for the configuration. If it is set to NA, a count of one is assumed for every
               row. This allows untabulated data to be processed. \code{cnts} must be a vector.
              }

  \item{runs}{Number of samples to be drawn.}
  \item{sig.item}{Indicator of significance in the result table (sig.z,sig.chisq,sig.perli,sig.zl,
sig.zl.corr). Do not forget to set the proper parameters for the CFA if sig.perli,sig.zl or
sig.zl.corr are to be used!}
  \item{\dots}{Parameters to be to relayed to the CFA}
}
\details{
 Takes 'runs' samples and does as many CFAs while counting how many times this configuration 
was  considered to be significant. 

Repeated-measures CFAs (mcfa) are not provided.

This is a heuristic method rather than a strict test of significance since there 
is no adjustment for multiple testing whatsoever. The advantage  is a more reliable picture
compared to splitting the original data, doing a CFA, and checking if the configurations re-appear 
in a CFA with the other half of the data. 
}
\value{
  \item{cnt.antitype}{Number of antiypes}
  \item{cnt.type }{Number of types}
  \item{pct.types}{Number of types in percent}
  \item{cnt.sig}{Number of significant results}
  \item{pct.cnt.sig}{Number of significant results in percent}
}
\references{ 
   Lautsch, E., von Weber S. (1995) Methoden und Anwendungen der Konfigurationsfrequenzanalyse 
   Psychologie und Medizin, Beltz Psychologie Verlagsunion
}
\author{Stefan Funke <s.funke@t-online.de>}

\note{\code{bcfa()} performs many CFAs which are by themselves slow, so the execution can
be \bold{very} time-consuming, especially if a sufficiently high value for \code{runs} was selected}

\seealso{ \code{\link{cfa}}, \code{\link{scfa}}}



\examples{
# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts<-trunc(runif(250)*10)
bcfa(configs,counts,runs=25) 
}

\keyword{multivariate}
\keyword{htest}

