#' Initialise functional splitting object elements
#' 
#' @param cl vector of clusters
#' @param top_set_genes set of genes modules used
#' @param best_set_genes_go details about the genes used
#' @param cliqueRes 
#' 
#' @return list of functional_split object elements init
#' @author Celine Barlier
init_obj <- function(cl,top_set_genes,best_set_genes_go,cliqueRes){
    #Clustering
    dfClust <- data.frame("Cell"=names(cl),"Cluster"=cl,"pathString"=rep(NA,length(cl)),"finalOutput"=rep(NA,length(cl)))
    dfClust$pathString[which(dfClust$Cluster == 1)] <- "allCells/1"
    dfClust$pathString[which(dfClust$Cluster == 2)] <- "allCells/2"
    dfClust$finalOutput[which(dfClust$Cluster == 1)] <- "0"
    dfClust$finalOutput[which(dfClust$Cluster == 2)] <- "1"
    #Set of genes for each split
    setGenesHC <- list()
    setGenesHC[["1|0"]] <- top_set_genes
    #GO enrichment for each split
    goHC <- list()
    goHC[["1|0"]] <- best_set_genes_go$GO[1]
    #CliqueTarget info
    cliqueTargets <- list()
    cliqueTargets[["1|0"]] <- cliqueRes[[best_set_genes_go$numCliques[1]]]
    
    return(list(dfClust,setGenesHC,goHC,cliqueTargets))
  }