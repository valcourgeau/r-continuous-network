
pckgs_to_load <- c(
  "igraph",
  "Rcpp",
  "stats",
  "prophet",
  "statnet",
  "intergraph",
  "visNetwork",
  "magrittr",
  "mlVAR",
  "graphicalVAR",
  "magrittr",
  "RTransferEntropy",
  "EntropyEstimation",
  "CompQuadForm",
  "JumpTest",
  "parallel",
  "RColorBrewer",
  "colorspace",
  "pracma"
)

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

check.packages(pckgs_to_load)
