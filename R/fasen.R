CHECK_EXAMPLE <- FALSE

FasenRegression <- function(data){
  data <- as.matrix(data)
  data_without_first <- data[-1,]
  data_without_last <- data[-nrow(data),]
  sub <- t(data_without_last) %*% data_without_last
  return(t(data_without_first) %*% data_without_last %*% solve(sub))
}

if(CHECK_EXAMPLE){
  FasenRegression(CleanData(data = df_load)$remainders)
}