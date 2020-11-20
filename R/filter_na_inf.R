#' Filter NA and Inf values
#' @author Alexey Samosyuk
filter_na_inf <-
  function(mtx)
  {
    colsS = log1p(Matrix::colSums(mtx))
    rowsS = log1p(Matrix::rowSums(mtx))
    
    goodCols = !is.infinite(colsS)&!is.na(colsS)&(colsS>0)
    goodRows = !is.infinite(rowsS)&!is.na(rowsS)&(rowsS>0)
    
    percC = (table(goodCols)/ncol(mtx))["TRUE"]
    percR = (table(goodRows)/nrow(mtx))["TRUE"]
    
    if(is.na(percR))
    {
      mtx = mtx[,goodCols]
    }else if(is.na(percC))
    {
      mtx = mtx[goodRows,]
    }else if(percC>percR)
    {
      mtx = mtx[,goodCols]
    }else{
      mtx = mtx[goodRows,]
    }
    
    
    return(mtx)
  }