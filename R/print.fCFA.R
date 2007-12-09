#print method for fCFA
print.fCFA <- function(x,...)
{
  devVec <- x$dev.val
  chisqVec <- x$chisq.val
  chidfVec <- x$df.val
  pvalueVec <- x$p.val
  strucMat <- x$struc.mat

  final.table <- data.frame(as.vector(round(devVec,2)),as.vector(round(chisqVec,2)),
  as.vector(round(chidfVec,2)),as.vector(round(pvalueVec,5)))
  dimnames(final.table)[[2]] <- c("LR","X^2","df","p")
  dimnames(final.table)[[1]] <- paste("Step",0:(dim(final.table)[1]-1))
  cat("Results of fCFA-fit: \n")
  cat("\n")
  print(final.table)
  cat("\n")
  cn <- colnames(strucMat)
  svec <- apply(strucMat,1, function(x) cn[x==1])
  ex <- data.frame(svec,x$typevec)
  dimnames(ex)[[2]] <- c("Excluded Cell","Type/Antitype")
  print(ex)
  cat("\n")
}