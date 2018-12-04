

genRdmSymmetricGraphs <- function(d, p.link){
  if(p.link <= 0.0 | p.link > 1.0)
    stop("p.link should be between 0 and 1.")
  
  
  bins <- sample(c(0,1), size = d*(d-1)/2, 
                 prob = c(1-p.link, p.link), replace = T)
  mat <- matrix(0, nrow = d, ncol = d)
  
  mat[upper.tri(mat, diag = F)] <- bins
  mat <- mat + t(mat)
  diag(mat) <- 1
  return(mat)
}

genRdmSymmetricGraphs(d = 8, p.link = 0.3)

genRdmAssymetricGraphs <- function(d, p.link){
  if(p.link <= 0.0 | p.link > 1.0)
    stop("p.link should be between 0 and 1.")
  
 
  return(matrix(
    sample(c(0,1), size = d*d, prob = c(1-p.link, p.link), replace = T),
    ncol=d))
}

genRdmAssymetricGraphs(d = 100, p.link = 0.15)

