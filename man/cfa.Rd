\name{cfa}
\alias{cfa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Analysis of configuration frequencies}
\description{
Determines the frequency of all combinations of variable values (their 
configurations) in comparison with their expected frequency calculated from 
the marginals and displays them in the order of decreasing chi-squared. In addition, a global chi squared is calculated.
}
\usage{
cfa(configmatrix, cntvector, 
    descending=TRUE, sort.on.chisq=TRUE, sort.on.n=FALSE,ignore.na=FALSE, 
    binom.test=FALSE, binom.test.limit=10, 
    bonferroni.p.z=T, bonferroni.alpha=0.05, 
    lehmacher=F, holm.alpha=0.01, verbose=FALSE)
}

\arguments{
 \item{configmatrix}{ Dataframe with the variables to be analyzed}
 \item{cntvector}{ Vector of counts (1 if the data are not aggregated)}
 \item{descending}{ Output in the order of decreasing chi squared}
 \item{sort.on.chisq}{ Sort output on chi squared}
 \item{sort.on.n}{ Sort output on the frequency of the configurations}
 \item{ignore.na}{ Ignore (casewise) missing data in the configurations }
 \item{binom.test}{ Perform binomial test for each configuration}
 \item{binom.test.limit}{ Maximum count (frequency) for which a binomial test 
is performed}
 \item{bonferroni.p.z}{ Bonferroni-adjust in sig. test for each configuration}
 \item{bonferroni.alpha}{ Alpha to be adjusted }
 \item{lehmacher}{ Perform Lehmacher's test with Holm's correction}
 \item{holm.alpha}{ Alpha to be adjusted according to Holm}
 \item{verbose}{ Long output}
}
\details{Each variable must have at least two different values and may have 
more (extension of the classical CFA). 

\code{configmatrix} must consist of  at least two variables (columns). 
Factors and numbers are both accepted (the numbers are internally converted 
to factors). \code{cntvector} must be numeric.

Counts should be at least = 5 for the chi squared test to be 
reliable but when using the CFA as a purely heuristic tool counts of 0 are 
possible. 

A z-approximation is used for the test of significance of 
configuration frequency.

If the data are not aggregated, i.e., there are several entries with the
same configuration (the same contents of a row in \code{configmatrix}) 
the counts of these configurations are added and one entry with the summed 
up count is generated, replacing the original entries.
}
\value{
  A list with class \code{"cfa"} contains the tabular results and the overall parameters
  \item{Row names}{Configuration}
  \item{n}{Frequency (count) of this configuration}
  \item{pct}{Relative Frequency (count) of this configuration}
  \item{expected}{Expected Frequency (count) of this configuration calculated from the marginals}
  \item{Q}{Coefficient of pronouncedness of the configuration. Varies between 0 and 1}
  \item{chisq}{Chi squared for the given configuration}
  \item{z}{z-approximation}
  \item{p}{p(z)}
  \item{sig(p(z))}{ 1: significant 0: not significant (limit is Bonferroni-adjusted by default}
  \item{Overall chi squared}{Overall chi squared for the entire table}
  \item{p(chi squared) }{p(chi squared) for the entire table }
  \item{Degrees of freedom }{Degrees of freedom for chi squared test of the entire table}
  \item{Total n}{Sum of all counts}
}
\references{
Krauth J., Lienert G. A. (1973, Reprint 1995)
Die Konfigurationsfrequenzanalyse (KFA) und ihre Anwendung in Psychologie und Medizin
Beltz Psychologie Verlagsunion

Eye, A. von (1990)
Introduction to configural frequency analysis. The search for types and anti-types in cross-classification. 
Cambride 1990
}
\author{ Stefan Funke <funke@attglobal.net>}
\note{The CFA is very useful a a heuristic tool when large numbers of categorical variables are to be screened because there is only one table of results instead of a multi-dimensional crosstabulation or a large amounts of sub-tables generated from it}

\section{WARNING }{The program is implemented in R itself rather than a compiled library and therefore slow. In most cases the input is a pre-aggregated table and speed is no problem because the \code{configmatrix} is small. There are no hard-coded limits in the program so even large tables can be processed but this will take time and memory. The outout table can be very wide if the levels of the factors variables are long 
strings so \code{options(width=..)} may need to be adjusted} 

\seealso{ \code{\link{mcfa}}, \code{\link{hier.cfa}}, \code{\link{boot.cfa}}
}

\examples{
library(cfa)
data(cfadat)
cfa(cfadat[c("gender","married","children")],cfadat["count"],verbose=T)
}
\keyword{misc}%-- one or more ...
