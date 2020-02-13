
CleanData <- function(data, frequency=24, s.window=24, t.window=24, ...){
  stl_cleaned <- lapply(1:ncol(data),
        FUN = function(i){
          stl(ts(data[,i], frequency = frequency), s.window = s.window, t.window = t.window, ...)
        }
  )
  names(stl_cleaned) <- colnames(data)
  remainders <- vapply(stl_cleaned, function(x){as.numeric(x$time.series[,3])}, rep(0, length(data[,1])))
  colnames(remainders) <- colnames(data)
  std.dev <- apply(remainders, MARGIN = 2, sd)
  # remainders <- t(apply(remainders, 1, function(x){x/std.dev}))
  return(list(stl_obj=stl_cleaned, remainders=remainders, std.dev=std.dev))
}