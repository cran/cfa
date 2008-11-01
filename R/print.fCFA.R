`print.fCFA` <-
function(x,...)
{
  cat("\n")
  cat("Results of fCFA-fit: \n")
  cat("\n")
  print(apply(x$restable, 2, round, 4))
  cat("\n")
  cn <- colnames(x$struc.mat)
  svec <- apply(x$struc.mat,1, function(xx) cn[xx==1])
  ex <- data.frame(svec,x$typevec)
  dimnames(ex)[[2]] <- c("Excluded Cell","Type/Antitype")
  print(ex)
  cat("\n")
}

