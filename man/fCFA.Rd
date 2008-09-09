\name{fCFA}
\alias{kvCFA}
\alias{fCFA}
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m.i}{Vector of observed frequencies.}
  \item{X}{Design Matrix of the base model.}
  \item{tabdim}{Vector of table dimensions.}
# \item{restype}{Residual type to be used; 
  \code{"stdPearson"} for standardized Pearson,
  \code{"Pearson"} for Pearson residuals.}
  \item{alpha}{Significance level.}
}

\value{
  \item{resstep}{Log-linear solution for each step: The first list element indicates the current iteration,
  the second the design matrix, the third the expected frequencies, the fourth is a vector composed of 
  the deviance, Chi-squared-value, df, and p-value.}
  \item{dev.val}{Deviance values for each step.}
  \item{chisq.val}{Chi-squared values for each step.}
  \item{df.val}{Degrees of freedom for each step.}
  \item{p.val}{P-values for each step.}
  \item{struMat}{Design vectors that blank out cells (structural part of the model).}
  
}
\references{

von Eye, A., and Mair, P. (2007). A functional approach to configural frequency analysis. 
Austrian Journal of Statistics, forthcoming.

Kieser, M., and Victor, N. (1999). Configural frequency analysis (CFA) revisited: A new look at an old approach. 
Biometrical Journal, 41, 967-983. 

}
\author{Patrick Mair, Alexander von Eye}

\examples{

#Functional CFA for a internet terminal usage data set by Wurzer (An application of configural frequency analysis: Evaluation of the
#usage of internet terminals, 2005, p.82)
dd <- data.frame(a1=gl(3,4),b1=gl(2,2,12),c1=gl(2,1,12))
X <- model.matrix(~a1+b1+c1,dd,contrasts=list(a1="contr.sum",b1="contr.sum",c1="contr.sum"))
ofreq <- c(121,13,44,37,158,69,100,79,24,0,26,3)
tabdim <- c(3,2,2)

res1 <- fCFA(ofreq, X, tabdim=tabdim)
res1
summary(res1)


#Kieser-Vector CFA for Children's temperament data from von Eye (Configural Frequency Analysis, 2002, p. 192) 
dd <- data.frame(a1=gl(3,9),b1=gl(3,3,27),c1=gl(3,1,27))
X <- model.matrix(~a1+b1+c1,dd,contrasts=list(a1="contr.sum",b1="contr.sum",c1="contr.sum"))
ofreq <- c(3,2,4,23,23,6,39,33,9,11,29,13,19,36,19,21,26,18,13,30,41,12,14,23,8,6,7)
tabdim <- c(3,3,3)

res2 <- kvCFA(ofreq, X, tabdim=tabdim)
res2
summary(res2)
}

\keyword{models}
