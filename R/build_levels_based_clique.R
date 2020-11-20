#' RECURSIVE FUNCTION: Build the child level: split in two the parent
#'
#' @param dfClust dataframe used to build the final object
#' @param scm single cell matrix normalized
#' @param scm_nm single cell matrix not normalized/filtered
#' @param firstLevel True if level 1, F if level > 1
#' @param clStop vector of indices to determine when to stop 
#' @param lev level number tracking
#' @param tfs vector of TFs names or IDs
#' @param setGenesHC list of set of genes used for each splitting
#' @param goHC list of GO for each set of genes used for each splitting
#' @param gda data.frame of GO slim annotations
#' @param cliqueTargets list of each antagonistic clique & targets used for the splitting
#' @param adjMethod multiple test correction to perform
#' @param cutoff p-adjusted p-value cutoff used for the enrichment
#'
#' @return list(dfClust,setGenesHC,goHC,cliqueTargets)
#' @author Celine Barlier
build_levels_based_clique <-
  function(dfClust,scm,scm_nm,firstLevel,clStop,lev=2,tfs,setGenesHC,goHC,gda,cliqueTargets,adjMethod,cutoff,qtarget,qprobInt,posRatio){
    
    #Check how many clusters - to know the level : level = # clusters / 2
    nc_actual <- length(unique(dfClust$Cluster))
    
    #For the level - each branch, splitted in two
    for (i in seq(1,nc_actual)) {
      if(!i %in% clStop){
        r <- build_branch_module(dfClust,i,(max(dfClust$Cluster) + 1),scm,scm_nm,clStop,tfs,setGenesHC,goHC,gda,cliqueTargets,adjMethod,cutoff,qtarget,qprobInt,posRatio)
        dfClust <- as.data.frame(r[[1]])
        clStop <- r[[2]]
        setGenesHC <- r[[3]]
        goHC <- r[[4]]
        cliqueTargets <- r[[5]]
      }else{
        #if i in clStop: add cell name to pathString if not already present
        if(length(which(str_detect(dfClust$pathString[which(dfClust$Cluster == i)],dfClust$Cell[which(dfClust$Cluster == i)]) == TRUE)) != length(dfClust$Cell[which(dfClust$Cluster == i)])){
          dfClust$pathString[which(dfClust$Cluster == i)] <- paste(dfClust$pathString[which(dfClust$Cluster == i)],dfClust$Cell[which(dfClust$Cluster == i)],sep = "/")
        }else{
          dfClust <- dfClust
        }
      }
    }
    
    #If all the cells are not found in the pathString > means that we can go deeper in the tree construction
    if(length(which(str_detect(dfClust$pathString,dfClust$Cell) == TRUE)) != length(dfClust$Cell)){
      return(build_levels_based_clique(dfClust,scm,scm_nm,F,clStop,lev = lev + 1,tfs,setGenesHC,goHC,gda,cliqueTargets,adjMethod,cutoff,qtarget,qprobInt,posRatio))
    }else{
      return(list(dfClust,setGenesHC,goHC,cliqueTargets))
    }
  }
