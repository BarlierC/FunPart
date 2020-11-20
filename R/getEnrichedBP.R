#' GO enrichment of a set of genes
#'
#' @param dataset single-cell dataset (cells in columns / genes in rows)
#' @param setGenes vector of genes (from one antagonistic module)
#' @param gda data.frame of GO slim annotations
#' @param adjMethod multiple test correction to perform
#' @param cutoff p-adjusted p-value cutoff used for the enrichment
#' 
#' @return data.frame of enriched BP 
#' @author Celine Barlier
getEnrichedBP <-
  function(dataset,setGenes,gda,adjMethod,cutoff){
    org <- ""
    if(str_detect(setGenes[1], "^[:upper:]+$")){
      #If uppercase genes > Human
      org <- "human"
    }else{
      #If not > Mouse
      org <- "mouse"
    }
    
    deg <- setGenes
    if(org == "human"){
      MB2gene=gda[, c("GOID", "Gene")]
    }else{
      MB2gene=gda[, c("GOID", "Gene")]
    }
    
    MB2name=gda[, c("GOID", "GOterm")]
    
    #universe = rownames(dataset)
    
    res <- enricher(deg, TERM2GENE=MB2gene, TERM2NAME=MB2name,pvalueCutoff = 0.05, pAdjustMethod = adjMethod, minGSSize = 5, maxGSSize = 2000,universe = rownames(dataset))
    
    if(length(res)>0){
      res <- as.data.frame(res@result) 
      if(length(rownames(res)) > 0){
        res <- res[which(res$p.adjust < cutoff),] #filter p-adj greater than 5%
        if(length(rownames(res)) > 0){
          res <- res[which(res$Count > 1),] #filter if only 1 gene match the category
        }else{
          return(data.frame()) 
        }
      }
      return(res)
    }else{
      return(data.frame()) 
    }
  }
