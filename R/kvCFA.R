`kvCFA` <-
function(m.i, X, tabdim, alpha=0.05)
#freq...vector of observedfrequencies
#X...full design matrix (with intercept)
{
  fit <- FALSE
  strucmat <- diag(1,length(m.i))                          #starting matrix for removing cells
  i <- -1
  chisqVec <- NULL
  chidfVec <- NULL
  pvalueVec <- NULL
  strucMat <- NULL
  devVec <- NULL
  resList <- NULL
  typevec <- NULL

  while ((fit==FALSE) || i==length(m.i)) {
    i <- i+1
    result <- glm.fit(X,m.i,family=poisson())              #model fit
    efreq <- result$fitted.values                            #expected frequencies
    chisq <- sum((m.i-efreq)^2/efreq)                      #chi-square value
    chisqVec <- cbind(chisqVec,chisq)
    chidf <- result$df.residual                              #degrees of freedom
    chidfVec <- cbind(chidfVec,chidf)
    devVec <- cbind(devVec,result$deviance)
    pvalue <- 1-pchisq(result$deviance,chidf)                          #pvalue
    pvalueVec <- cbind(pvalueVec,pvalue)
    fitvec <- c(chisq,result$deviance,result$df.residual,pvalue)
    names(fitvec) <- c("LR","X^2","df","p") 
    Xdes <- as.matrix(X)
    colnames(Xdes) <- NULL
    rownames(Xdes) <- NULL
    rL <- list(desmat = Xdes, exp.freq = efreq, fitvec = fitvec)
    resList <- c(resList,rL)

    if (pvalue < alpha) {                                    #add structural coding vector
      fit <- FALSE
      fitstruc <- apply(strucmat,2,function(x) {
                                    Xnew <- cbind(X,x)       #new design matrix with removed cells
                                    result <- glm.fit(Xnew,m.i,family=poisson())           #model fit
                                    efreq <- result$fitted.values                            #expected frequencies
                                    chisq <- sum((m.i-efreq)^2/efreq)                      #chi-square value
                                    list(chisq,x)
                                    })
      fitstruc1 <- unlist(fitstruc,recursive=FALSE)
      chivec <-  unlist(fitstruc1[seq(1,length(fitstruc1),by=2)])    #chi-squared values
      indChiMin <- (1:length(chivec))[min(chivec)==chivec]           #position with minimum chi-squared value
      strucvec <- fitstruc1[seq(2,length(fitstruc1),by=2)]           #list with structural design vectors
      minStruc <- strucvec[[indChiMin]]                              #structural vector to be added
      strucMat <- rbind(strucMat,minStruc)
      X <- cbind(X,minStruc)
      TFmat <- strucmat==minStruc
      posDel <- (1:dim(TFmat)[2])[colSums(TFmat)==dim(TFmat)[1]]
      strucmat <- strucmat[,-posDel]
      if (m.i[minStruc==1] > efreq[minStruc==1]) {
        typ = "type"
      } else {
        typ = "antitype"
      }
      typevec <- c(typevec,typ)
    } else {
      fit <- TRUE
    }
  }
  
 # final.table <- data.frame(as.vector(round(devVec,2)),as.vector(round(chisqVec,2)),as.vector(round(chidfVec,2)),as.vector(round(pvalueVec,5)))
 # dimnames(final.table)[[2]] <- c("L","X2","df","p")
  
  if (is.null(strucMat)) stop("Base model fits! No types/antitypes found!\n",call.=FALSE)
  
  rownames(strucMat) <- paste("Step",1:(dim(strucMat)[1]))     #matrix with excluded elements
  td <- length(tabdim)
  gridlist <- tapply(tabdim[td:1],1:td,function(x) 1:x)
  namesmat <- expand.grid(gridlist)[,td:1]
  strvec <- apply(namesmat,1,paste,collapse="")
  colnames(strucMat) <- strvec
  
  
  restable <- data.frame(as.vector(devVec), as.vector(chisqVec), as.vector(chidfVec), as.vector(pvalueVec)) 
  dimnames(restable)[[2]] <- c("LR","X2","df","p")
  dimnames(restable)[[1]] <- paste("Step",0:(dim(restable)[1]-1))
  names(typevec) <- paste("Step",1:length(typevec))
  
  desmat <- resList[(length(resList)-2):length(resList)][[1]]               #design matrix
  
  lcomp <- unique(names(resList))
  nsteps <- dim(restable)[1]
  steplab <- paste("step", 1:nsteps, sep = "")
  names(resList) <- as.vector(t(outer(steplab, lcomp, paste, sep = "")))
  
  result <- list(restable = restable, design.mat = desmat,  struc.mat = strucMat, 
                 typevec = typevec, resstep = resList)

  class(result) <- "kvCFA"
  result
}

