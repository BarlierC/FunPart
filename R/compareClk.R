#' Identify if a clique is unique
#' @param v vectors of TFs
#' @param clk liste of positive cliques
#' @return true or false 
#' @author Celine Barlier
compareClk <- 
  function(v,clk,co){
    if(length(co)>0){clk <- clk[-co]}
    n <- foreach(i=seq(1,length(clk)),.combine = "c") %dopar% {
      length(intersect(names(clk[[i]]),names(v)))/length(names(clk[[i]]))
    }

    #Remove 100% = same clique
    n <- n[n<1]
    
    #If 80% similar with other cliques = not unique
    if(length(which(n>=0.7))>0){
      return(FALSE) #not unique
    }else{
      return(TRUE) #unique 
    }
}