\name{plot.bcfa}
\alias{plot.bcfa}
\title{Plotting method for a bcfa object}
\description{Plots an object of the class \code{bcfa}}

\usage{
\method{plot}{bcfa}(x,...)
}

\arguments{
  \item{x}{An object of the class \code{bcfa} which is returned by the function 
  \code{boot.cfa()}}
  \item{...}{Any arguments to be given to \code{plot}}
}
\details{
  Plots the number of cases considered significant vs. the number of cases considered
  to be a type (n > expected). 

  This is in some way like other plots of quality versus quantity.

  Configurations can be identified by left-clicking on them until the
  right mouse button is pressed. The labels of the configurations selected will be displayed
  in the text window.
}
\value{
Returns a vector of the configurations selected with their name set to the labels
}

\references{None - plots have been rarely used with the CFA}
\author{Stefan Funke <s.funke@t-online.de>}

\note{This function is usually invoked plotting an object returned by 
\code{bcfa}}

\seealso{ \code{\link{bcfa}}}

\examples{
# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts<-trunc(runif(250)*10)
plot(bcfa(configs,counts,runs=25))
}

\keyword{multivariate}
\keyword{htest}

