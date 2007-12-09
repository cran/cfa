summary.fCFA <- function(object,...)
{
  devVec <- object$dev.val
  chisqVec <- object$chisq.val
  chidfVec <- object$df.val
  pvalueVec <- object$p.val
  strucMat <- object$struc.mat

  final.table <- data.frame(as.vector(round(devVec,2)),as.vector(round(chisqVec,2)),
  as.vector(round(chidfVec,2)),as.vector(round(pvalueVec,5)))
  dimnames(final.table)[[2]] <- c("LR","X^2","df","p")
  dimnames(final.table)[[1]] <- paste("Step",0:(dim(final.table)[1]-1))
  cat("Results of fCFA-fit: \n")
  cat("\n")
  print(final.table)
  cat("\n")

  cat("Final log-linear model: \n")                                 #print out last model
  final <- object$resstep[(length(object$resstep)-2):length(object$resstep)]
  cat("Design Matrix:\n")                                           #design matrix
  colnames(final[[1]]) <- NULL
  print(final[[1]])
  cat("\n Expected frequencies:\n")
  print(round(final[[2]],2))
  cat("\n")
}