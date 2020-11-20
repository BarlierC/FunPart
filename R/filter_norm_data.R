#' Clean (and normalize) the data 
#' @param mtx single cell matrix (cells in columns / genes in rows)
#' @param norm perform Seurat's normalization (True) or not (False)
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
    
    #Keep genes expressed at least in 10% of the cells - remove the one expressed in more than 90% of the cells
    m <- as.data.frame(m)
    percent_expressed=((rowSums(m != 0))*100)/ncol(m)
    m=cbind(m,percent_expressed)
    m=m[m$percent_expressed>percExp,]
    m$percent_expressed<-NULL
    m <- as.matrix(m)
    
    #Remove genes lowly expressed - to get macrophages res
    exp=(rowSums(m != 0))
    lowExp=exp[exp<quantile(exp,qExp)] 
    m <- m[which(!rownames(m) %in% names(lowExp)),]
    
    return(m)
  }