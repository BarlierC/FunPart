#' CALLED BY build_levels_based_clique() to build cluster for the two branches of the level x
#'
#' @param dfClust dataframe as Cell | Cluster
#' @param actuC number of actual cluster
#' @param newC integer used to define the next cluster
#' @param scm single cell RNA-seq matrix
#' @param scm_nm single cell matrix not normalized/filtered
#' @param clStop vector of indices to determine when to stop
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
build_branch_module <-
  function(dfClust,actuC,newC,scm,scm_nm,clStop,tfs,setGenesHC,goHC,gda,cliqueTargets,adjMethod,cutoff,qtarget,qprobInt,posRatio){
    
    #If actuC not in clStop
    if(!actuC %in% clStop){
      
      #New ip
      scm <- as.data.frame(scm) 
      scm_nm <- as.data.frame(scm_nm)
      nip <- scm[,which(colnames(scm) %in% dfClust$Cell[which(dfClust$Cluster == actuC)])]
      nip_nm <- scm_nm[,which(colnames(scm_nm) %in% dfClust$Cell[which(dfClust$Cluster == actuC)])]
      
      #Keep genes expressed in more than 10%
      percent_expressed=((rowSums(nip != 0))*100)/ncol(nip)
      nip=cbind(nip,percent_expressed)
      nip=nip[nip$percent_expressed>10,]
      nip$percent_expressed<-NULL

      #If at least 10 cells
      if(length(colnames(nip)) > 10){
        
        #Build network & get set of genes based on cliques
        cliqueRes <- find_mutual_cliques(as.matrix(nip),tfs,qtarget) 
        
        #If we do not find clique - stop
        if(length(cliqueRes) == 0){
          clStop <- c(clStop,actuC)
          return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
        }else{
          #clique selection  - Perform GO enrichment 
          tmp <- data.frame("SetGenes"=rep("",length(cliqueRes)))
          for (i in seq(1,length(cliqueRes))){tmp$SetGenes[i] <- paste(c(paste(unique(c(names(cliqueRes[[i]][[1]]),as.character(c(unlist(cliqueRes[[i]][[1]]))))),collapse=",")),c(paste(unique(c(names(cliqueRes[[i]][[2]]),as.character(c(unlist(cliqueRes[[i]][[2]]))))),collapse=",")),sep="&")} #Structure: gene_set1,gene_set1,gene_set1&gene_set2,gene_set2,gene_set2
          best_set_genes_go <- runGOanalysis_listSetGenes(as.data.frame(nip_nm),tmp,gda,adjMethod,cutoff)
          
          #If NA enrichment - stop
          if(length(best_set_genes_go$FC)==0){
            clStop <- c(clStop,actuC)
            return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
          }else{
            #Top set of genes
            top_set_genes <- c(strsplit(strsplit(best_set_genes_go$SetGenes[1],split="&")[[1]][1],",")[[1]],strsplit(strsplit(best_set_genes_go$SetGenes[1],split="&")[[1]][2],",")[[1]])
            rm(tmp)
            #Heatmap
            mhc <- nip[which(rownames(nip) %in% top_set_genes),]
            
            #If no SD == 0, use it
            if(length(which(colSums(mhc) == 0))==0){
              use <- TRUE
              clNum <- 1
            }else if(length(which(colSums(mhc) == 0))>0 & length(best_set_genes_go$SetGenes) > 1){
              #If SD == 0 but other cliques exists, test them
              rc <- test_new_clique_from_list(scm,best_set_genes_go,incS=2,cellNum=5)
              use <- rc[[1]]
              mhc <- rc[[2]]
              top_set_genes <- rc[[3]]
              clNum <- rc[[4]]
              cl <- rc[[5]]
            }else{
              #If SD == 0 and no other clique - stop
              use <- FALSE
            }

              if(use){
                #Clustering
                h <- hierarchical_clust(mhc)
                if(length(h)>1){
                  #Get clusters
                  h <- as.hclust(h)
                  cl <- cutree(h,k=2)
                  if(length(unique(cl)) == 1){
                    #If only one cluster > we stop
                    clStop <- c(clStop,actuC)
                    return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
                  }else if(table(cl)[1] < 5 || table(cl)[2] < 5 & length(best_set_genes_go$SetGenes) > 1){
                    #If less than 5 cells but other cliques exists - try another one
                    rc <- test_new_clique_from_list(nip,best_set_genes_go,incS = 2, 5)
                    usethen <- rc[[1]]
                    mhc <- rc[[2]]
                    top_set_genes <- rc[[3]]
                    cl <- rc[[5]]
                    if(usethen){
                      cl[which(cl == 2)] <- newC
                      cl[which(cl == 1)] <- actuC
                      dfClust$Cluster[which(dfClust$Cell %in% names(cl))] <- cl
                      dfClust$pathString[which(dfClust$Cluster == newC)] <- paste(dfClust$pathString[which(dfClust$Cluster == newC)],'2',sep = "/")
                      dfClust$pathString[which(dfClust$Cluster == actuC)] <- paste(dfClust$pathString[which(dfClust$Cluster == actuC)],'1',sep = "/")
                      
                      setGenesHC[[paste(dfClust$finalOutput[which(dfClust$Cluster == newC)][1],"_1|0",sep="")]] <- top_set_genes
                      
                      goHC[[paste(dfClust$finalOutput[which(dfClust$Cluster == newC)[1]],"_1|0",sep="")]] <- best_set_genes_go$GO
                      
                      cliqueTargets[[paste(dfClust$finalOutput[which(dfClust$Cluster == newC)[1]],"_1|0",sep="")]] <- cliqueRes[[best_set_genes_go$numCliques[clNum]]]
                      
                      dfClust$finalOutput[which(dfClust$Cluster == newC)] <- paste(dfClust$finalOutput[which(dfClust$Cluster == newC)],'1',sep = "_")
                      dfClust$finalOutput[which(dfClust$Cluster == actuC)] <- paste(dfClust$finalOutput[which(dfClust$Cluster == actuC)],'0',sep = "_")
                      return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
                    }else{
                      #we stop for this branch
                      clStop <- c(clStop,actuC)
                      return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
                    }
                  }else if(table(cl)[1] >= 5 & table(cl)[2] >= 5){
                    cl[which(cl == 2)] <- newC
                    cl[which(cl == 1)] <- actuC
                    dfClust$Cluster[which(dfClust$Cell %in% names(cl))] <- cl
                    dfClust$pathString[which(dfClust$Cluster == newC)] <- paste(dfClust$pathString[which(dfClust$Cluster == newC)],'2',sep = "/")
                    dfClust$pathString[which(dfClust$Cluster == actuC)] <- paste(dfClust$pathString[which(dfClust$Cluster == actuC)],'1',sep = "/")
                    
                    setGenesHC[[paste(dfClust$finalOutput[which(dfClust$Cluster == newC)][1],"_1|0",sep="")]] <- top_set_genes
                    
                    goHC[[paste(dfClust$finalOutput[which(dfClust$Cluster == newC)[1]],"_1|0",sep="")]] <- best_set_genes_go$GO
                    
                    cliqueTargets[[paste(dfClust$finalOutput[which(dfClust$Cluster == newC)[1]],"_1|0",sep="")]] <- cliqueRes[[best_set_genes_go$numCliques[clNum]]]
                    
                    dfClust$finalOutput[which(dfClust$Cluster == newC)] <- paste(dfClust$finalOutput[which(dfClust$Cluster == newC)],'1',sep = "_")
                    dfClust$finalOutput[which(dfClust$Cluster == actuC)] <- paste(dfClust$finalOutput[which(dfClust$Cluster == actuC)],'0',sep = "_")
                    return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
                  }else{
                    #we stop for this branch
                    clStop <- c(clStop,actuC)
                    return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
                  }  
                }else{
                  #we stop for this branch
                  clStop <- c(clStop,actuC)
                  return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
                }
              }else{
                #we stop for this branch
                clStop <- c(clStop,actuC)
                return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
              }
            }
        }
      }else{
        #Less than 10 cells to analyze - stop 
        clStop <- c(clStop,actuC)
        return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
      }   
    }
    return(list(dfClust,clStop,setGenesHC,goHC,cliqueTargets))
  }
