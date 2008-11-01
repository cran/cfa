\name{plot.mcfa}
\alias{plot.mcfa}
\title{ Plotting method for a mcfa object}
\description{
Plots an object of the class \code{mcfa}
}
\usage{
\method{plot}{mcfa}(x,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{An object of the class \code{mcfa} which is returned by the function 
  \code{cfa()} (rather than \code{mcfa()}) which performs a repeated measures CFA (two or more columns of counts)}
  \item{...}{Any arguments to be used by \code{plot}}

}
\details{
Plots chi squared vs. the sum of all counts for this configuration which indicates 
pronouncedness of the configuration vs. practical importance. 
Configurations can be identified by left-clicking on them until the
right mouse button is pressed. The labels of the configurations selected will be displayed
in the text window.
}
\value{Returns a list of the labels of the configurations selected.
}

\references{None - plots have been rarely used with the CFA}
\author{Stefan Funke <s.funke@t-online.de>}

\note{This function is usually invoked plotting an object returned by 
\code{cfa}}

\seealso{ \code{\link{cfa}}, \code{\link{mcfa}}}


\examples{
# Some random configurations:
configs<-cbind(c("A","B")[rbinom(250,1,0.3)+1],c("C","D")[rbinom(250,1,0.1)+1],
          c("E","F")[rbinom(250,1,0.3)+1],c("G","H")[rbinom(250,1,0.1)+1])
counts1<-trunc(runif(250)*10)
counts2<-trunc(runif(250)*10)

plot(cfa(configs,cbind(counts1,counts2)))

}
\keyword{multivariate}
\keyword{htest}
