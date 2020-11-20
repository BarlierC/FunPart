#' Set names
#' @author Alexey Samosyuk
set_names <-
  function(mtx)
  {
    if(is.null(colnames(mtx)))
    {
      colnames(mtx) = 1:ncol(mtx)
    }
    
    if(is.null(rownames(mtx)))
    {
      rownames(mtx) = 1:nrow(mtx)
    }
    
    mtx = mtx[rownames(mtx) %in% names(which(table(rownames(mtx))==1)),]
    
    return(mtx)
  }