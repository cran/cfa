\name{plot.hcfa}
\alias{plot.hcfa}
\title{Plotting method for a hcfa object}
\description{
Plots an object of the class \code{hcfa}
}
\usage{
\method{plot}{hcfa}(x,...)
}
\arguments{
  \item{x}{An object of the class hcfa}
  \item{...}{Any arguments to be used by \code{plot}}
}
\details{ A dotchart is generated which plots chi squared vs. the order of the configuration 
(i.e. the number of configuration variables it contains).
}
\value{
Returns NULL.
}

\references{None - plots have been rarely used with the CFA}
\author{Stefan Funke <s.funke@t-online.de>}

\note{This function is usually invoked plotting an object returned by 
\code{hcfa}}

\seealso{ \code{\link{cfa}}, \code{\link{hcfa}}}


\examples{
#configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
#          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
#counts<-trunc(runif(250)*10)
#plot(hcfa(configs,counts))
}

\keyword{multivariate}
\keyword{htest}
