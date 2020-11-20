#' Clean the data matrix
#' @author Alexey Samosyuk
filter_mtx <-
  function(mtx, th = 1)
  {
    mtx = as(mtx, "dgCMatrix")
    mtx[mtx<0]=0
    mtx = filter_na_inf(mtx)
    
    mtx = set_names(mtx)
    
    return(mtx)
  }