\name{cfadat}
\alias{cfadat}
\non_function{}
\title{ Simple cfa test data}
\usage{data(cfadat)}
\description{
The \code{cfadat} data frame has 8 rows and 4 columns.
There are 3 columns containing the configurations and 1 column containing 
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
\item{count}{
a numeric vector containing the counts
}
}
}
\details{
A dataset for cfa must contain 2 or more columns with factors or numbers 
indication the configuration and one column indicating the counts
}
\source{
  The data were freely invented to illustrate the function of the cfa
}
\examples{
data(cfadat)
cfa(cfadat[1:3],cfadat[4])
}
\keyword{datasets}
