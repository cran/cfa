\name{lcfa}
\alias{lcfa}
\title{ Test according to Lindner}
\description{
  Performs an hypergeometrical, nonparametrical exact test of significance according to Lindner (1984)
}
\usage{
lcfa(m, n, Nt, k, P.XisM0=TRUE, P.XatMostM0=TRUE, P.XatLeastM0=TRUE)
}

\arguments{ Usually, just one of the p values available is required for hypothesis testing. It may save time to exclude the rest.
 \item{m}{ Observed frequency of the configuration tested}
 \item{n}{ Marginal sums of the parameters realized in the configuration to be tested (vector) }
 \item{Nt}{ Sample size of configurations (N total) }
 \item{k}{ Number of parameters }
 \item{P.XisM0}{ Return p(X=m). }
 \item{P.XatMostM0}{ Return p(X<=m). }
 \item{P.XatLeastM0}{ Return p(X>=m). }
}


\value{ Returns p values for the hypergeometrical test according to Linder (1984)
 \item{pHXequalsM0}{ p(X=m). The probability of obtaining exactly m occurrences of the configuration given n and Nt }
 \item{pHXsmallerEqualM0}{ p(X<=m). The probability of obtaining at most m occurrences of the configuration given n and Nt }
 \item{pHlargerEqualM0}{ p(X>=m).The probability of obtaining at least m occurrences of the configuration given n and Nt }
 \item{timed.required}{ The computation time spent (milliseconds). }
}

\references{
Lindner, K.: Eine exakte Auswertungsmethode zur Konfigurationsfrequenzanalyse [An exact procedure for
the configural frequency analysis]. Psycholog Beitraege 26, 393?415 (1984)

Harloff, Joachim, An efficient algorithm for Lindners test
(configural frequency analysis), Qual Quant DOI 10.1007/s11135-011-9499-9
     }

\author{J. Harloff <joachimharloff@joachimharloff.de> }

\note{ The test according to Lindner is boosted by the package multicore if available.
     }



\seealso{ \code{\link{cfa}}}

\examples{
lk<-4 # number of parameters
ln<-c(59,57,59,58) # marginal sums of the parameters realized in the configuration to be tested
lNt<-116 # sample size of configurations
lm0<-16 # observed frequency of the configuration tested
lcfa(lm0,ln,lNt,lk) # return all p values.


}

\keyword{multivariate}
\keyword{htest}
