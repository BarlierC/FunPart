#' Build the functional splitting object
#'
#' @param t list of results (cliques,GO,clusters)
#' @param scm single cell matrix 
#' 
#' @return functional_split object
#' @author Celine Barlier
build_object <- function(t,scm){
    #If no cliques is found
    if(length(t) == 0){
      #Build empty final object
      clust <- rep(1,length(colnames(scm)))
      names(clust) <- colnames(scm)
      
      setG <- list()
      setG[["Clust1"]] <- list()
      
      functE <- list()
      functE[["Clust1"]] <- list()
      
      cliqT <- list()
      cliqT[["Clust1"]] <- list()
      
      #Object
      functional_res <- list("data"=scm,"clust"=clust,"genesets"=setG,"function"=functE,"cliques"=cliqT)
      class(functional_res) <- "functionalSplitting"
    }else{
      #Build final object with results
      clust <- c(t[[1]]$finalOutput)
      names(clust) <- t[[1]]$Cell
      
      setG <- list()
      setG[["Clust1"]] <- t[[2]]
      
      functE <- list()
      functE[["Clust1"]] <- t[[3]]
      
      cliqT <- list()
      cliqT[["Clust1"]] <- t[[4]]
      
      #Object
      functional_res <- list("data"=scm,"clust"=clust,"genesets"=setG,"function"=functE,"cliques"=cliqT)
      class(functional_res) <- "functionalSplitting"
    }
    return(functional_res)
  }