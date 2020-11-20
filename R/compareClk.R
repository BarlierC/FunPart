#' Identify if a clique is unique
#' @param v vectors of TFs
#' @param clk liste of positive cliques
#' @return true or false 
#' @author Celine Barlier
compareClk <- 
  function(v,clk){
    
    n <- sapply(clk,function(x){
      length(intersect(names(x),names(v)))/length(names(x))
    })
    
    #Remove 100% = same clique
    n <- n[n<1]
    
    #If 80% similar with other cliques = not unique
    if(length(which(n>=0.7))>0){
      return(FALSE) #not unique
    }else{
      return(TRUE) #unique 
    }
  }