#' Update FunPart object structure
#' Add modules & functionalenrich slots to provide more details
#' @param FunPartObj FunPart object
#' @param m Initial single cell matrix (cells in columns, genes in rows)
#' @param gda Enrichment ressource file 
#' @param padjMeth P-adjusted correction method
#' @param padj p-adjusted value cutoff
#' @return FunPart object updated slots
#' @author Celine Barlier
fillEnrichFunPartObj <-
function(FunPartObj,m,gda,padjMeth=adjMethod,padj=cutoff){
  if(length(FunPartObj$cliques$Clust1)>0){
    #Add New slots enrichment
    FunPartObj$functionalenrich <- list()
    FunPartObj$modules <- list()
    levName <- c(names(FunPartObj$cliques$Clust1))
    #For each module
    for(j in seq(1,(length(FunPartObj$cliques$Clust1)))){
      m1 <- paste(getEnrichedBPManual(m,unique(c(names(FunPartObj$cliques$Clust1[[j]]$C1),unique(unlist(FunPartObj$cliques$Clust1[[j]]$C1)))),gda,adjMethod,cutoff),collapse = ",")
      m2 <- paste(getEnrichedBPManual(m,unique(c(names(FunPartObj$cliques$Clust1[[j]]$C2),unique(unlist(FunPartObj$cliques$Clust1[[j]]$C2)))),gda,adjMethod,cutoff),collapse = ",")
      #Put in the object
      FunPartObj$modules[[levName[j]]]$M1 <- paste(unique(c(names(FunPartObj$cliques$Clust1[[j]]$C1),unique(unlist(FunPartObj$cliques$Clust1[[j]]$C1)))),collapse = ",")
      FunPartObj$modules[[levName[j]]]$M2 <- paste(unique(c(names(FunPartObj$cliques$Clust1[[j]]$C2),unique(unlist(FunPartObj$cliques$Clust1[[j]]$C2)))),collapse = ",")
      FunPartObj$functionalenrich[[levName[j]]]$M1 <- m1
      FunPartObj$functionalenrich[[levName[j]]]$M2 <- m2
    }
    #Remove old Enrichment slot
    FunPartObj$`function` <- NULL
    return(FunPartObj)
  }else{
    print("No functional state identified in this object")
    FunPartObj$functionalenrich <- list()
    FunPartObj$modules <- list()
    FunPartObj$`function` <- NULL
    return(FunPartObj)
  }
}