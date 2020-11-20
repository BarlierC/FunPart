#' Hierarchical clustering using Pheatmap package
#'
#' @param mhc single cell RNA-seq with genes being the ones of the best antagonistic and functionally enriched modules
#'
#' @return dendrogram object
#' @author Celine Barlier
hierarchical_clust <-
  function(mhc){
    heatmap <- pheatmap(log2(mhc+1),clustering_distance_cols="correlation") 
    h <- as.dendrogram(as.hclust(heatmap$tree_col))
    return(h)
  }