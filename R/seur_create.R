#' Create a Seurat object
#' @author Alexey Samosyuk
seur_create <-
  function(mtx)
  {
    object = CreateSeuratObject(counts = mtx, names.delim = "no_delimiters_plz!")
    object[["orig.ident"]] = get_types(colnames(mtx))
    object[["curr.ident"]] = object@active.ident
    return(object)
  }