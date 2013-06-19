\name{PXisM}
\alias{PXisM}
\title{ Test according to Lindner}
\description{
  Performs a test of significance according to Lindner
}
\usage{
PXisM(m,n,Nt,k)
}

\arguments{
 \item{m}{Observed frequency of the observation tested}
 \item{n}{Marginal sums of the parameters realized in the configuration to be tested (vector) }
 \item{Nt}{Sample size of configurations }
 \item{k}{Number of parameters}

}


\value{ returns p for the test according to Linder

}

\references{
Lindner, K.: Eine exakte Auswertungsmethode zur Konfigurationsfrequenzanalyse [An exact procedure for
the configural frequency analysis]. Psycholog Beitraege 26, 393?415 (1984)

Harloff, Joachim, An efficient algorithm for Lindners test
(configural frequency analysis), Qual Quant DOI 10.1007/s11135-011-9499-9
     }

\author{J. Harloff <oachimharloff@joachimharloff.de> }

\note{ The test according to Lindner requires the packages parallel. All other parts of cfa do not.
     }



\seealso{ \code{\link{cfa}}}

\examples{
# Does not work with windows since there is no parallel for it
if (require(parallel)) {
lk<-4 # number of parameters
ln<-c(59,57,59,58) # marginal sums of the parameters realized in the configuration to be tested
lNt<-116 # sample size of configurations
lm0<-16 # observed frequency of the configuration tested


# New algorithm
starttime=proc.time()
pHXsmallerequalM0<-sum(unlist(mclapply(0:lm0,PXisM,ln,lNt,lk)))
pHXequalM0<-PXisM(lm0,ln,lNt,lk)
pHlargerequalM0<-sum(unlist(mclapply(lm0: min(ln),PXisM,ln,lNt,lk)))
stoptime<-proc.time()
list(pHXsmallerequalM0=pHXsmallerequalM0,pHXequalM0=pHXequalM0,pHlargerequalM0=pHlargerequalM0,
timed.required=stoptime-starttime)

# End of the new algorithm
}

}

\keyword{multivariate}
\keyword{htest}
