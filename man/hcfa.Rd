\name{hcfa}
\alias{hcfa}
\title{Hierachical analysis of configuration frequencies}
\description{Recursively eliminates one variable in the configuration to generate all possible
sub-tables and performs a global chi-squared-test on them
}
\usage{
hcfa(configs, cnts)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{configs}{ Contains the configurations. This can be a dataframe or a matrix. The dataframe can contain numbers,
              characters, factors or booleans. The matrix can consist of numbers, characters or booleans 
              (factors are implicitely re-converted to numerical levels). There must be >=3 columns.
            }
  \item{cnts}{ Contains the counts for the configuration. If it is set to NA, a count of one is assumed for every
               row. This allows untabulated data to be processed. \code{cnts} can be a vector or a matrix/dataframe with >=2 columns.
              }
}
\details{The hierarchical CFA assists in the selection of configuration variables by showing
which variables contribute the most to the variability. If eliminating a variable does not
markedly decrease the global chi squared the variable is likely to be redundant, provided there are
no extraneous reasons for retaining it. 

The output is in decreasing order of chi squared so the most useful combinations of 
variables come first.
  
}
\value{
  \item{chisq}{Global chi squared}
  \item{df}{Degrees of freedom for this subtable}
  \item{order}{Order (number of configuration variables)}
}
\references{Lautsch, E., von Weber S. (1995) Methoden und Anwendungen der Konfigurationsfrequenzanalyse 
     in Psychologie und Medizin, Beltz Psychologie Verlagsunion}
\author{Stefan Funke <s.funke@t-online.de>}

\note{The p for the test of significance ist provided by the print method}

\seealso{\code{\link{cfa}}, \code{\link{scfa}}, \code{\link{mcfa}}}

\examples{# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],
c("C","D")[rbinom(250,1,0.1)+1],
c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts<-trunc(runif(250)*10)
hcfa(configs,counts) 
}

\keyword{multivariate}
\keyword{htest}
