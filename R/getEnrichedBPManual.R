#' Get all the BPs enriched (used to fill in the new structure of FunPart object)
#' @param dataset Initial single cell matrix (cells in columns, genes in rows)
#' @param setGenes Enrichment ressource file 
#' @param gda Enrichment ressource file
#' @param padjMeth P-adjusted correction method
#' @param padj p-adjusted value cutoff
#' @return vector of functional enrichment
#' @author Celine Barlier
getEnrichedBPManual <- 
function(dataset,setGenes,gda,padjMeth=adjMethod,padj=cutoff){
  dataset <- dataset[which(!str_detect(rownames(dataset),"^ERCC")),]
  dataset <- dataset[which(!str_detect(rownames(dataset),"^mt")),]
  dataset <- as.matrix(dataset)
  deg <- setGenes
  MB2gene=gda[, c("GOID", "Gene")]
  MB2name=gda[, c("GOID", "GOterm")]
  res <- clusterProfiler::enricher(deg, TERM2GENE=MB2gene, TERM2NAME=MB2name,pvalueCutoff = 0.05, pAdjustMethod = padjMeth, minGSSize = 5, maxGSSize = 2000,universe = rownames(dataset))
  res <- as.data.frame(res@result)
  res <- res[which(res$p.adjust < padj),]
  res <- res[which(res$Count > 1),]
  resbp <- paste(res$ID,res$Description,sep=":")
  return(resbp)
}