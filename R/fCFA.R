`fCFA` <-
function(m.i, X, tabdim, alpha = 0.05)
{

#tabdim... a vector with table dimensions, e.g., tabdim = c(3,2,3)

  #if (alpha.cor == "Bonferroni") alpha <- alpha/length(m.i)

  resList <- NULL
  chisqVec <- NULL
  chidfVec <- NULL
  pvalueVec <- NULL
  strucMat <- NULL
  devVec <- NULL
  typevec <- NULL
  i <- -1

  fit <- FALSE
  startdim <- dim(X)[2]
  res <- "Pearson"
  #bonf <- length(m.i)
  while (fit==FALSE)
  {
    i <- i+1
    result <- glm.fit(X,m.i,family=poisson())              #model fit
    efreq <- result$fitted.values                            #expected frequencies

    if (res == "Pearson") {
      stres <- (m.i-efreq)/sqrt(efreq)                     #standardized Pearson residuals
      stresa <- abs(stres)                                   #standardized residuals (absolute)
      pvec <- 1-pnorm(stresa)
    }
    ##if (res == "Poisson") {                                  #standardized Poisson residuals
##      w.i <- result$weights
##      W <- diag(w.i)
##      Hat <- (mroot(W))%*%X%*%(solve(t(X)%*%W%*%X))%*%t(X)%*%(mroot(W))
##      h.i <- diag(Hat)
##      stres <- (m.i-efreq)/sqrt(efreq*(1-h.i))
##      stresa <- abs(stres)
##      pvec <- 1-pnorm(stresa)
##    }
##    if (res == "Deviance") {
##      d.i <- 2*m.i*log(m.i/efreq)
##      stres <- sqrt(abs(d.i))*sign(m.i-efreq)
##      stresa <- abs(stres)
##      pvec <- 1-pnorm(stresa)
##    }

    chisq <- sum((m.i-efreq)^2/efreq)                      #chi-square value
    chisqVec <- cbind(chisqVec,chisq)
    devVec <- cbind(devVec,result$deviance)
    chidfVec <- cbind(chidfVec,result$df.residual)
    pvalue <- 1-pchisq(result$deviance,result$df.residual)                          #pvalue
    pvalueVec <- cbind(pvalueVec,pvalue)
    fitvec <- c(chisq,result$deviance,result$df.residual,pvalue)
    names(fitvec) <- c("LR","X^2","df","p") 
    Xdes <- as.matrix(X)
    colnames(Xdes) <- NULL
    rownames(Xdes) <- NULL
    rL <- list(desmat = Xdes, exp.freq = efreq, fitvec = fitvec)
    resList <- c(resList,rL)
    
    if (pvalue < alpha) {
      strucvec <- as.integer(stresa==max(stresa))            #cells to be blanked out
      strucMat <- rbind(strucMat,strucvec)
      X <- cbind(X,strucvec)                                 #new design matrix
      fit <- FALSE
      if (m.i[strucvec==1] > efreq[strucvec==1]) {
        typ = "type"
      } else {
        typ = "antitype"
      }
      typevec <- c(typevec,typ)
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

 restable <- data.frame(as.vector(devVec), as.vector(chisqVec), as.vector(chidfVec), as.vector(pvalueVec)) 
  dimnames(restable)[[2]] <- c("LR","X^2","df","p")
  dimnames(restable)[[1]] <- paste("Step",0:(dim(restable)[1]-1))
  names(typevec) <- paste("Step",1:length(typevec))
  
  desmat <- resList[(length(resList)-2):length(resList)][[1]]               #design matrix
  
  lcomp <- unique(names(resList))
  nsteps <- dim(restable)[1]
  steplab <- paste("step", 1:nsteps, sep = "")
  names(resList) <- as.vector(t(outer(steplab, lcomp, paste, sep = "")))
  
  result <- list(restable = restable, design.mat = desmat,  struc.mat = strucMat, 
                 typevec = typevec, resstep = resList)
 
  class(result) <- "fCFA"
  
  result
}

