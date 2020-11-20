#' Get types
#' @author Alexey Samosyuk
get_types <-
  function(arr, sep = ".")
  {
    return(factor(sapply(arr, function(cell) { return(unlist(strsplit(cell, sep, fixed = T))[1])})))
  }