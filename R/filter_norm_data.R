#' Clean (and normalize) the data 
#' @param mtx single cell matrix (cells in columns / genes in rows)
#' @param norm perform Seurat's normalization (True) or not (False)
#' @param percExp percentage of cells in which the gene has to be expressed to be considered
#' @param qExp quantile above which the genes are considered (QC step)
#' 
#' @return cleaned (and normalized) matrix
#' @author Celine Barlier
filter_norm_data <-
  function(mtx,norm,percExp,qExp){
    mtx = Matrix(as.matrix(mtx),sparse = T)
    mtx = filter_mtx(mtx)
    
    obj = seur_create(mtx)
    rm(mtx)
    
    #Normalize data if norm = T
    if(norm){
      obj = NormalizeData(obj) 
    }
    
    #matrix of normalized data
    m = obj@assays$RNA@data
    
    #QC cells: Remove outlier cells if any 
    n <- m
    n[n>0] <- 1
    x <- colSums(n)
    names(x) <- colnames(n)
    #Hampel filter for outliers
    threshold <- median(x) - 3 * mad(x) 
    #Statistical test on the potential outliers
    if(length(which(x<threshold))>0){
      stattest <- rosnerTest(x,k = length(which(x<threshold)))
      #Get outliers if any
      if(length(stattest$all.stats$Outlier)>0){
        out <- names(x[which(x %in% stattest$all.stats$Value[which(stattest$all.stats$Outlier == T)])])
      }
      if(length(out) > 0){
        #Remove outliers
        m <- m[,which(!colnames(m) %in% out)]
      }
    }
    
    #Keep genes expressed at least in x% of the cells 
    m <- as.data.frame(m)
    percent_expressed=((rowSums(m != 0))*100)/ncol(m)
    m=cbind(m,percent_expressed)
    m=m[m$percent_expressed>percExp,]
    m$percent_expressed<-NULL
    m <- as.matrix(m)
    
    #QC genes: Remove genes too lowly expressed compared to the dataset
    exp=(rowSums(m != 0))
    lowExp=exp[exp<quantile(exp,qExp)] 
    m <- m[which(!rownames(m) %in% names(lowExp)),]
    
    return(m)
  }