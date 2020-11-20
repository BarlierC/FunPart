#' Test different cliques to be used if one get drop out all along the criteria filters
#' 
#' @param cl vector of clusters
#' @param top_set_genes set of genes modules used
#' @param best_set_genes_go details about the genes used
#' 
#' @return list of functional_split object elements init
#' @author Celine Barlier
test_new_clique_from_list <- function(scm,best_set_genes_go,incS=2,cellNum=5){
  checkSD <- FALSE
  use <- FALSE
  
  while(checkSD == F){
    top_set_genes <- c(strsplit(strsplit(best_set_genes_go$SetGenes[incS],split="&")[[1]][1],",")[[1]],strsplit(strsplit(best_set_genes_go$SetGenes[incS],split="&")[[1]][2],",")[[1]])
    mhc <- scm[which(rownames(scm) %in% top_set_genes),]
    #If no SD = 0
    if(length(which(colSums(mhc) == 0)) == 0){
      #If the splitting lead to clusters with more than 5 cells
      h <- hierarchical_clust(mhc)
      h <- as.hclust(h)
      cl <- cutree(h,k=2)
      #If one of the two clusters contains less than cellNum cells, try another gene module
      if(table(cl)[1] >= cellNum & table(cl)[2] >= cellNum){
        checkSD <- TRUE #break the boucle
        use <- TRUE #use the module found
      }else if(incS == length(best_set_genes_go$SetGenes)){
        checkSD <- TRUE #break the boucle
        use <- FALSE #no clique corresponds
      }
    }else if(incS == length(best_set_genes_go$SetGenes)){
      checkSD <- TRUE #break the boucle
      use <- FALSE #no clique corresponds
    }
    incS <- incS + 1
  }
  
  return(list(use,mhc,top_set_genes,incS))
}