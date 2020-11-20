#' Run BP enrichment for a list of set of genes
#' 
#' @param dataset single-cell dataset (cells in columns / genes in rows)
#' @param listSetGenes list of identified antagonistic modules
#' @param gda data.frame of GO slim annotations
#' @param adjMethod multiple test correction to perform
#' @param cutoff p-adjusted p-value cutoff used for the enrichment
#'
#' @return df of cliques + GO + FC enrichment
#' @export
#' @author Celine Barlier
runGOanalysis_listSetGenes <-
  function(dataset,listSetGenes,gda,adjMethod,cutoff){
    
    enrichment_setGenes_1 <- foreach(a=1:length(rownames(listSetGenes)),.export=c("runGOanalysis","getEnrichedBP","listSetGenes","gda","adjMethod","cutoff"),.packages = c("stringr","clusterProfiler","Seurat")) %dopar% {
      runGOanalysis(dataset,strsplit(strsplit(listSetGenes$SetGenes[a],split="&")[[1]][1],",")[[1]],gda,adjMethod,cutoff)
    }
    
    enrichment_setGenes_2 <- foreach(b=1:length(rownames(listSetGenes)),.export=c("runGOanalysis","getEnrichedBP","listSetGenes","gda","adjMethod","cutoff"),.packages = c("stringr","clusterProfiler","Seurat")) %dopar% {
      runGOanalysis(dataset,strsplit(strsplit(listSetGenes$SetGenes[b],split="&")[[1]][2],",")[[1]],gda,adjMethod,cutoff)
    }
    
    listSetGenes$GO <- rep(NA,length(listSetGenes$SetGenes))
    listSetGenes$FC <- rep(NA,length(listSetGenes$SetGenes))
    listSetGenes$numCliques <- seq(1,length(listSetGenes$SetGenes))
    
    #For each antagonistic cliques get the GO ID + term + FC enrichment score
    for (i in seq(1,length(rownames(listSetGenes)))) {
      r_set1 <- as.data.frame(enrichment_setGenes_1[[i]])
      r_set2 <- as.data.frame(enrichment_setGenes_2[[i]])
      #If the two antagonistic set of genes are enriched
      if(length(rownames(r_set1)) > 0 & length(rownames(r_set2)) > 0){
        
        #Antagonistic set gene 1
        r_set1 <- r_set1[order(r_set1$Count,decreasing = T),]
        
        #Antagonistic set gene 2
        r_set2 <- r_set2[order(r_set2$Count,decreasing = T),]
        
        listSetGenes$GO[i] <- paste(paste(r_set1$ID[1],r_set1$Description[1],sep=" "),paste(r_set2$ID[1],r_set2$Description[1],sep=" "),sep="&")
        listSetGenes$FC[i] <- (sum(as.numeric(sapply(r_set1$GeneRatio, function(x) eval(parse(text=x)))))/length(r_set1$ID)) * (sum(as.numeric(sapply(r_set2$GeneRatio, function(x) eval(parse(text=x)))))/length(r_set2$ID))
      }
    }
    
    #Order by FC enrichment
    listSetGenes <- listSetGenes[order(listSetGenes$FC,decreasing = T),] 
    
    #Remove the ones to NA
    listSetGenes <- listSetGenes[!is.na(listSetGenes$FC),]
    
    #Return the top clique
    if(!is.na(listSetGenes$FC[1])){
      return(listSetGenes)
    }else{
      return(data.frame())
    }
  }