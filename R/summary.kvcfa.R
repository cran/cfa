`summary.kvCFA` <-
function(object,...)
{
  cat("Results of KV-CFA-fit: \n")
  cat("\n")
  print(apply(object$restable, 2, round, 4))
  cat("\n")

  cat("Final log-linear model: \n")                                 #print out last model
  final <- object$resstep[(length(object$resstep)-2):length(object$resstep)]
  cat("Design Matrix:\n")                                           #design matrix
  print(object$design.mat)
  cat("\n Expected frequencies:\n")
  print(round(final[[2]],2))
  cat("\n")
}

