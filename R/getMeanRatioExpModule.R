#' Get mean ratio cell expressing TFs of the module
#' 
#' @param m normalized/cleaned single cell expression matrix
#' @param moduletfs vector of TFs names belonging to the same module
#' @param cells cells to use to calculate the expression ratio
#' 
#' @return mean ratio of cells expressing moduletfs
#' @author Celine Barlier
getMeanRatioExpModule <- function(m,moduletfs,cells){
  m <- as.matrix(m)
  m <- m[which(rownames(m) %in% moduletfs),which(colnames(m) %in% cells)]
  mb <- m
  #Binary matrix: expressed - 1, not expressed - 0
  mb[mb>0] <- 1
  mb[mb<0] <- 0
  r <- rowSums(mb)/ncol(mb) #Compute ratio
  return(mean(r))
}
