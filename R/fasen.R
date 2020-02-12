
FasenRegression <- function(data){
  data <- as.matrix(data)
  data_without_first <- data[-1,]
  data_without_last <- data[-nrow(data),]
  t(data_without_last)
  sub <- t(data_without_last) %*% data_without_last
  
  return(t(data_without_first) %*% data_without_last %*% solve(sub))
}

FasenRegression(df_load)
