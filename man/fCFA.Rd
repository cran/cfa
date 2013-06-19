\name{fCFA}

\alias{fCFA}
\alias{kvCFA}
\alias{print.fCFA}
\alias{print.kvCFA}
\alias{summary.fCFA}
\alias{summary.kvCFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stepwise CFA approaches}
\description{These CFA methods detect and eliminate stepwise types/antitypes
cells by specifying an appropriate contrast in the design matrix.
The procedures stop when model fit is achieved. Functional CFA (fCFA) uses a residual
criterion, Kieser-Victor CFA (kvCFA) a LR-criterion.
}
\usage{
fCFA(m.i,  X, tabdim, alpha = 0.05)
kvCFA(m.i, X, tabdim, alpha = 0.05)
}


\arguments{
  \item{m.i}{Vector of observed frequencies.}
  \item{X}{Design Matrix of the base model.}
  \item{tabdim}{Vector of table dimensions.}
  \item{alpha}{Significance level.}
}

\value{
  \item{restable}{Fit results for each step}
  \item{design.mat}{Final design matrix}
  \item{struc.mat}{Structural part of the design matrix for each step}
  \item{typevec}{Type or antitype for each step}
  \item{resstep}{Design matrix, expected frequency vector, and fit results for each step}  
}
\references{

von Eye, A., and Mair, P. (2008). A functional approach to configural frequency analysis. 
Austrian Journal of Statistics, 37, 161-173.

Kieser, M., and Victor, N. (1999). Configural frequency analysis (CFA) revisited: A new look at an old approach. 
Biometrical Journal, 41, 967-983. 

}
\author{Patrick Mair, Alexander von Eye}

\examples{

#Functional CFA for a internet terminal usage data set by Wurzer 
#(An application of configural frequency analysis: Evaluation of the
#usage of internet terminals, 2005, p.82)
dd <- data.frame(a1=gl(3,4),b1=gl(2,2,12),c1=gl(2,1,12))
X <- model.matrix(~a1+b1+c1,dd,contrasts=list(a1="contr.sum",b1="contr.sum",
    c1="contr.sum"))
ofreq <- c(121,13,44,37,158,69,100,79,24,0,26,3)
tabdim <- c(3,2,2)

res1 <- fCFA(ofreq, X, tabdim=tabdim)
res1
summary(res1)


# Kieser-Vector CFA for Children's temperament data from 
# von Eye  (Configural Frequency Analysis, 2002, p. 192) 
dd <- data.frame(a1=gl(3,9),b1=gl(3,3,27),c1=gl(3,1,27))
X <- model.matrix(~a1+b1+c1,dd,contrasts=list(a1="contr.sum",
    b1="contr.sum",c1="contr.sum"))
ofreq <- c(3,2,4,23,23,6,39,33,9,11,29,13,19,36,19,21,26,18,13,30,
         41,12,14,23,8,6,7)
tabdim <- c(3,3,3)

res2 <- kvCFA(ofreq, X, tabdim=tabdim)
res2
summary(res2)
}

\keyword{models}
