
genLevyCovMatrix <- function(d){
  A <- matrix(runif(d^2)*2-1, ncol=d) 
  Sigma <- t(A) %*% A
}

genRdmSymmetricGraphs <- function(d, p.link, theta_1 = 1, theta_2 = 1){
  if(p.link <= 0.0 | p.link > 1.0)
    stop("p.link should be between 0 and 1.")
  
  
  bins <- sample(c(0,1), size = d*(d-1)/2, 
                 prob = c(1-p.link, p.link), replace = T)
  mat <- matrix(0, nrow = d, ncol = d)
  
  mat[upper.tri(mat, diag = F)] <- bins
  mat <- theta_1*(mat + t(mat))
  diag(mat) <- theta_2
  return(mat)
}

genRdmSymmetricGraphs(d = 8, p.link = 0.3, theta_1 = 1, theta_2 = 1)

genRdmAssymetricGraphs <- function(d, p.link, theta_1 = 1, theta_2 = 1){
  if(p.link <= 0.0 | p.link > 1.0)
    stop("p.link should be between 0 and 1.")
  
  res <- theta_1 * matrix(
    sample(c(0,1), size = d*d, prob = c(1-p.link, p.link), replace = T),
    ncol=d)
  diag(res) <- theta_2
  return(res)
}

genRdmAssymetricGraphs(d = 100, p.link = 0.15)

genRdmAssymmetricExpBinGraphs <- function(d, p.link, rate){
  # By Dr. Luitgard Veraart
  L1<-matrix(rexp(n=d*d, rate=rate), nrow=d, ncol=d)
  L2<-matrix(rbinom(n=d*d, size=1, p=p.link), nrow=d, ncol=d)
  L<-round(L1*L2, 0)
  diag(L)<-0
  rownames(L)<-1:d
  colnames(L)<-1:d
  return(L)
}

genRdmAssymmetricExpBinGraphs(d = 5, p.link = 0.25, rate = 0.1)

