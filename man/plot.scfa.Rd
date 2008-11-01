\name{plot.scfa}
\alias{plot.scfa}
\title{ Plotting method for a scfa object}
\description{Plots an object of the class \code{scfa}
}
\usage{
\method{plot}{scfa}(x,...)
}
\arguments{
  \item{x}{An object of the class \code{scfa} which is returned by the function 
  \code{cfa()} (rather than \code{scfa()}) which performs a simple CFA (one column of counts)}
  \item{...}{Any arguments to be used by \code{plot}}
}
\details{
Plots chi squared vs. n which indicates pronouncedness of the configuration vs. 
practical importance. Configurations can be identified by left-clicking on them until the
right mouse button is pressed. The labels of the configurations selected will be displayed
in the text window.
}
\value{
  Returns a list of the labels of the configurations selected.
}

\references{None - plots have been rarely used with the CFA}
\author{Stefan Funke <s.funke@t-online.de>}

\note{This function is usually invoked plotting an object returned by 
\code{cfa}}

\seealso{ \code{\link{cfa}}, \code{\link{scfa}}}


\examples{# library(cfa) if not yet loaded
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts<-trunc(runif(250)*10)
plot(cfa(configs,counts))
}

\keyword{multivariate}
\keyword{htest}


