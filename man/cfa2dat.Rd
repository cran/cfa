\name{cfa2dat}
\alias{cfa2dat}
\non_function{}
\title{Example data for multiple cfa}
\usage{data(cfa2dat)}
\description{
The \code{cfa2dat} data frame has 8 rows and 5 columns.
There are 3 columns containing the configurations and 2 columns containing 
the counts
}
\format{
This data frame contains the following columns:
\describe{
\item{gender}{
a factor with levels
\code{f} 
\code{m} 
}
\item{married}{
a factor with levels
\code{mno} 
\code{myes} 
}
\item{children}{
a factor with levels
\code{cno} 
\code{cyes} 
}
\item{count1950}{
a numeric vector with the first set of counts
}
\item{count1990}{
a numeric vector with the second set of counts
}
}
}
\details{
A dataset for mcfa must contain 2 or more columns with factors or numbers 
indicating the configuration and two columns indicating the counts
}
\source{
The data were freely invented for illustrative purposes only
}
\examples{
data(cfa2dat)
mcfa(cfa2dat[1:3],cfa2dat[4:5],verbose=T)
}
\keyword{datasets}
