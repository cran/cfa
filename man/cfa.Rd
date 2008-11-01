\name{cfa}
\alias{cfa}
\title{ Analysis of configuration frequencies}
\description{
  This is the main function which will call scfa() und mcfa() as required to handle the simple and the multiple cfa.
}
\usage{
cfa(cfg, cnts=NA, sorton="chisq", sort.descending=TRUE, format.labels=TRUE, 
    casewise.delete.empty=TRUE, 
    binom.test=FALSE, exact.binom.test=FALSE, exact.binom.limit=10, 
    perli.correct=FALSE, lehmacher=FALSE, lehmacher.corr=TRUE, 
    alpha=0.05, bonferroni=TRUE)
}
\arguments{
  \item{cfg}{ Contains the configurations. This can be a dataframe or a matrix. The dataframe can contain numbers,
              characters, factors, or booleans. The matrix can consist of numbers, characters, or booleans 
              (factors are implicitely re-converted to numerical levels). There must be >=3 columns.
            }
  \item{cnts}{ Contains the counts for the configuration. If it is set to NA, a count of one is assumed for every
               row. This allows untabulated data to be processed. \code{cnts} can be a vector or a matrix/dataframe with >=2 columns.
              }
  \item{sorton}{ Determines the sorting order of the output table. Can be set to \code{chisq}, \code{n}, or \code{label}.}
  \item{sort.descending}{ Sort in descending order}
  \item{format.labels}{ Format the labels of the configuration. This makes to output wider but it will increase the
        readability.}
  \item{casewise.delete.empty}{ If set to TRUE all configurations containing a NA in any column will be deleted. Otherwise
        NA is handled as the string "NA" and will appear as a valid configuration.}
  \item{binom.test}{ Use z approximation for binomial test.}
  \item{exact.binom.test}{ Do an exact binomial test.}
  \item{exact.binom.limit}{ Maximum n for which an exact binomial test is performed (n >10 causes p to become inexact).}
  \item{perli.correct}{ Use Perli's correction for multiple test.}
  \item{lehmacher}{ Use Lehmacher's correction for multiple test.}
  \item{lehmacher.corr}{ Use a continuity correction for Lehmacher's correction.}
  \item{alpha}{Alpha level}
  \item{bonferroni}{Do Bonferroni adjustment for multiple test (irrelevant for Perli's and Lehmacher's test).}
}
\details{
  The cfa is used to sift large tables of nominal data. Usually it is used for dichotomous 
  variables but can be extended to three or more possible values. There should be at least three configuration variables
  in \code{cfg} - otherwise a simple contigency table would do.
  All tests of significance are two-sided: They test for both types or antitypes, i.e. if n is significantly larger or
  smaller than the expected value.
  The usual caveats for testing contigency tables apply. If a configuration has a n <5 an exact test
  should be used. As an alternative the least interesting configuration variable can be left out (if
  it is not essential) which will automatically increase the n for the remaining configurations.
}
\value{ Some of these elements will only be returned when the corresponding argument
        in the function call has been set. The relation is obvious due to corresponding names.
  \item{table}{The cfa output table}
  \item{table["label"]}{Label for the given configuration}
  \item{table["n"]}{Observed n for this configuration} 
  \item{table["expected"]}{Expected n for this configuration} 
  \item{table["Q"]}{Coefficient of pronouncedness (varies between 0 and 1)} 
  \item{table["chisq"]}{Chi squared for the given configuration} 
  \item{table["p.chisq"]}{p for the chi squared test} 
  \item{table["sig.chisq"]}{Is it significant (will Bonferroni-adjust if argument \code{bonferroni} 
  is set)} 
  \item{table["z"]}{z-approximation for chi squared}
  \item{table["p.z"]}{p of z-test}
  \item{table["sig.z"]}{Is it significant (will Bonferroni-adjust if argument \code{bonferroni} 
  is set)?}
  \item{table["x.perli"]}{Statistic for Perli's test}
  \item{table["sig.perli"]}{Is it significant (this is designed as a multiple test)?}
  \item{table["zl"]}{z for Lehmacher's test}
  \item{table["sig.zl"]}{Is it significant (this is designed as a multiple test)?}  
  \item{table["zl.corr"]}{z for Lehmacher's test (with continuity correction)}
  \item{table["sig.zl.corr"]}{Is it significant (this is designed as a multiple test)?}
  \item{table["p.exact.bin"]}{p for exact binomial test} 
  \item{summary.stats}{Summary stats for entire table}
  \item{summary.stats["totalchisq"]}{Total chi squared} 
  \item{summary.stats["df"]}{Degrees of freedom}
  \item{summary.stats["p"]}{p for the chi squared test}
  \item{summary.stats["sum of counts"]}{Sum of all counts}
  \item{levels}{Levels for each configuration. Should all be 2 for the bivariate case} 
}
\references{ Krauth J., Lienert G. A. (1973, Reprint 1995) Die
     Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in
     Psychologie und Medizin. Beltz Psychologie Verlagsunion
     
     Lautsch, E., von Weber S. (1995) Methoden und Anwendungen der Konfigurationsfrequenzanalyse in
     Psychologie und Medizin. Beltz Psychologie Verlagsunion

     Eye, A. von (1990) Introduction to configural frequency analysis.
     The search for types and anti-types in cross-classification. 
     Cambride 1990
     }
\author{Stefan Funke <s.funke@t-online.de>}
\note{There are no hard-coded limits in the program so even large tables can be processed. 
     The output table can be very wide if the levels of factors variables are long  strings so `options(width=..)'
     may need to be adjusted. 

     The object returned has the class scfa if a one-sample CFA was performed or 
     the class mcfa if a repeated-measures CFA was performed. \code{cfa()} decides which
     one is appropriate by looking at \code{cnts}: If it is a vector, it will do a
     simple CFA. If it is a dataframe or matrix with 2 or more columns, a repeated-measures
     CFA ist done.
     }

\section{WARNING }{Note than spurious "significant" configurations are likely to appear in very large tables.
     The results should therefore be replicated before they are accepted as real. \code{boot.cfa} can be helpful to check
     the results.} 

\seealso{ \code{\link{scfa}}, \code{\link{mcfa}}}

\examples{
# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts<-trunc(runif(250)*10)
cfa(configs,counts) 
}
\keyword{multivariate}
\keyword{htest}







