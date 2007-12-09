#library("mgcv")

#Functional CFA
fCFA <- function(n.i, X, tabdim, restype = "stdPearson", alpha = 0.05)
{

#tabdim... a vector with table dimensions, e.g., tabdim = c(3,2,3)
#n.i ... vector with observed frequencies
#restype ... type of residual ("stdPearson", "Pearson", "Deviance")
#X...design matrix with intercept

  ok <- match.arg(restype, c("stdPearson","Pearson"))
  #if (alpha.cor == "Bonferroni") alpha <- alpha/length(m.i)

  resList <- NULL
  chisqVec <- NULL
  chidfVec <- NULL
  pvalueVec <- NULL
  strucMat <- NULL
  devVec <- NULL
  i <- -1

  fit <- FALSE
  startdim <- dim(X)[2]
  typevec <- NULL
 
  while (fit==FALSE)
  {
    i <- i+1
    result <- glm.fit(X,n.i,family=poisson())              #model fit
    #result <- glm(m.i~X, family = poisson())  
    efreq <- result$fitted.values                            #expected frequencies

    if (restype == "Pearson") {
      stres <- (n.i-efreq)/sqrt(efreq)                     #standardized Pearson residuals
      stresa <- abs(stres)                                   #standardized residuals (absolute)
     # pvec <- 2*(1-pnorm(stresa))
    }
    if (restype == "stdPearson") {                                  #standardized Poisson residuals
      #design matrix noch den intercept weg!!
      w.i <- result$weights
      W <- diag(w.i)
      Hat <- (sqrt(W))%*%X%*%(solve(t(X)%*%W%*%X))%*%t(X)%*%(sqrt(W))
      h.i <- diag(Hat)
      stres <- rep(0,length(h.i))
      ind <- which(round(h.i,5) != 1)                     #ind: where ofreq != efreq
      stres[ind] <- (n.i[ind]-efreq[ind])/sqrt(efreq[ind]*(1-h.i[ind]))  #stand. Pearson residual
      stresa <- abs(stres)
     # pvec <- 2*(1-pnorm(stresa))
    }
    if (restype == "Deviance") {
      d.i <- 2*n.i*log(n.i/efreq)
      stres <- sqrt(abs(d.i))*sign(n.i-efreq)
      stresa <- abs(stres)
     # pvec <- 2*(1-pnorm(stresa))
    }

    chisq <- sum((n.i-efreq)^2/efreq)                      #chi-square value
    chisqVec <- cbind(chisqVec,chisq)
    devVec <- cbind(devVec,result$deviance)
    chidfVec <- cbind(chidfVec,result$df.residual)
    pvalue <- 1-pchisq(result$deviance,result$df.residual)                          #pvalue
    pvalueVec <- cbind(pvalueVec,pvalue)

    rL <- list(paste("Step ",i),X,efreq,c(result$deviance,chisq,result$df.residual,pvalue))
    resList <- c(resList,rL)
    if (pvalue < alpha) {
      strucvec <- as.integer((stresa==max(stresa)))            #cells to be blanked out
      print(efreq)
      if (n.i[strucvec==1] > efreq[strucvec==1]) {
        typ = "type"
      } else {
        typ = "antitype"
      }
      typevec <- c(typevec,typ) 
      strucMat <- rbind(strucMat,strucvec)
      X <- cbind(X,strucvec)                                 #new design matrix
      fit <- FALSE
    } else {
      fit <- TRUE
    }
    if (result$df.residual==0) {
      fit <- TRUE                   #no more df left
      warning("No more df left. Iteration is abandoned.")
    }
  }
  
  if (is.null(strucMat)) stop("Base model fits! No types/antitypes found!\n",call.=FALSE)
    
  rownames(strucMat) <- paste("Step",1:(dim(strucMat)[1]))     #matrix with excluded elements
  td <- length(tabdim)
  gridlist <- tapply(tabdim[td:1],1:td,function(x) 1:x)
  namesmat <- expand.grid(gridlist)[,td:1]
  strvec <- apply(namesmat,1,paste,collapse="")
  colnames(strucMat) <- strvec

  result <- list(resstep = resList, dev.val = devVec, chisq.val = chisqVec,
  df.val = chidfVec, p.val = pvalueVec, struc.mat = strucMat, typevec = typevec)
  
  class(result) <- "fCFA"
  
  result
}

#result:
#resstep ... results for each fCFA-step
#dev.val ... deviance values for each step
#chisq.val ... Chi-square values for each step
#df.val ... corresponding df
#p.val... p-values for each step