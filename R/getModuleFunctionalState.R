#' Process FunPart objects content into a table with all results
#' 
#' @param funpartobj FunPart output object
#' 
#' @return dataframe with the formatted results: modules, branch, type, tfs, genes and enrichment
#' @export
#' @author Celine Barlier
getModuleFunctionalState <- function(funpartobj){
  #Df
  dfres <- data.frame("Module"=NA,"Branch"=NA,"Type"=NA,"TFs"=NA,"Genes"=NA,"Enrichment"=NA)
  #Get functional states identified
  csnb <- unique(funpartobj$clust)
  #If splitting happened
  if(length(csnb)>1){
    #For each module level (split)
    lev <- names(funpartobj$cliques$Clust1)
    
    for (i in seq(1,length(lev)))
    {
      if(lev[i] == "1|0"){
        #First level
        mo <- c("1","0") #Branch modules to test
      }else{
        #Level second and more
        tmp <- strsplit(lev[i],split="_")[[1]]
        tmp <- tmp[-c(length(tmp))]
        tmp <- paste(tmp,collapse = "_")
        mo <- c(
          paste(tmp,"1",sep="_"),
          paste(tmp,"0",sep="_")
        )
      }

      #If branch is a leaf = functional state
      typeMod <- list()
      if(mo[1] %in% csnb & !mo[2] %in% csnb){
        #Branch 1 is a leaf = functional state and not branch 0
        cells_branch1 <- names(funpartobj$clust[which(funpartobj$clust == mo[1])])
        typeMod[["branch1"]] <- "Direct module"
        
        group_branch0 <- csnb[which(str_detect(csnb,paste("^",mo[2],sep="")) == TRUE)] #All cluster under branch 0
        cells_branch0 <- names(funpartobj$clust[which(funpartobj$clust %in% group_branch0)])
        typeMod[["branch0"]] <- "Intermediate module"
        
      }else if(mo[2] %in% csnb & !mo[1] %in% csnb){
        #Branch 0 is a leaf = functional state and not branch1
        cells_branch0 <- names(funpartobj$clust[which(funpartobj$clust == mo[2])])
        typeMod[["branch0"]] <- "Direct module"
        
        group_branch1 <- csnb[which(str_detect(csnb,paste("^",mo[1],sep="")) == TRUE)] 
        cells_branch1 <- names(funpartobj$clust[which(funpartobj$clust %in% group_branch1)])
        typeMod[["branch1"]] <- "Intermediate module"
        
      }else if(mo[1] %in% csnb & mo[2] %in% csnb){
        #Both are direct modules
        cells_branch1 <- names(funpartobj$clust[which(funpartobj$clust == mo[1])])
        typeMod[["branch1"]] <- "Direct module"
        
        cells_branch0 <- names(funpartobj$clust[which(funpartobj$clust == mo[2])])
        typeMod[["branch0"]] <- "Direct module"
        
      }else{
        group_branch1 <- csnb[which(str_detect(csnb,paste("^",mo[1],sep="")) == TRUE)] #All cluster under branch 1
        group_branch0 <- csnb[which(str_detect(csnb,paste("^",mo[2],sep="")) == TRUE)] #All cluster under branch 0
        #No branch is a leaf
        cells_branch1 <- names(funpartobj$clust[which(funpartobj$clust %in% group_branch1)])
        cells_branch0 <- names(funpartobj$clust[which(funpartobj$clust %in% group_branch0)])
        typeMod[["branch0"]] <- "Intermediate module"
        typeMod[["branch1"]] <- "Intermediate module"
      }
      
      print(i)
      
      #MODULE 1
      #Test M1 - belongs to branch 0
      m1_branch0 <- getMeanRatioExpModule(funpartobj$data,
                                          names(funpartobj$cliques$Clust1[[lev[i]]]$C1),
                                          cells_branch0)
      #Test M1 - belongs to branch 1
      m1_branch1 <- getMeanRatioExpModule(funpartobj$data,
                                          names(funpartobj$cliques$Clust1[[lev[i]]]$C1),
                                          cells_branch1)
  
      #MODULE 2
      #Test M2 - belongs to branch 0
      m2_branch0 <- getMeanRatioExpModule(funpartobj$data,
                                          names(funpartobj$cliques$Clust1[[lev[i]]]$C2),
                                          cells_branch0)
      #Test M2 - belongs to branch 1
      m2_branch1 <- getMeanRatioExpModule(funpartobj$data,
                                          names(funpartobj$cliques$Clust1[[lev[i]]]$C2),
                                          cells_branch1)
      
      #Comparisons
      if(m1_branch0 > m1_branch1 & m2_branch0 < m2_branch1){
        #M1 belongs to branch 0 & M2 belongs to branch 1
        #Branch 0
        dfres <- rbind(dfres,
                       data.frame("Module"="M1",
                                  "Branch"=mo[2],
                                  "Type"=typeMod[["branch0"]],
                                  "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C1),collapse = ","),
                                  "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C1)),collapse = ","),
                                  "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M1))
        #Branch 1
        dfres <- rbind(dfres,
                       data.frame("Module"="M2",
                                  "Branch"=mo[1],
                                  "Type"=typeMod[["branch1"]],
                                  "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C2),collapse = ","),
                                  "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C2)),collapse = ","),
                                  "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M2))
        
      }else if(m1_branch1 > m1_branch0 & m2_branch1 < m2_branch0){
        #M1 belongs to branch 1 & M2 belongs to branch 0
        #Branch 0
        dfres <- rbind(dfres,
                       data.frame("Module"="M1",
                                  "Branch"=mo[1],
                                  "Type"=typeMod[["branch1"]],
                                  "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C1),collapse = ","),
                                  "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C1)),collapse = ","),
                                  "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M1))
        #Branch 1
        dfres <- rbind(dfres,
                       data.frame("Module"="M2",
                                  "Branch"=mo[2],
                                  "Type"=typeMod[["branch0"]],
                                  "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C2),collapse = ","),
                                  "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C2)),collapse = ","),
                                  "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M2))
        
      }else{
        #Ambuiguity may arise for higher level intermediate levels, 
        #in this case check the higher cell exp ratio & attribute the modules accordingly
        rb <- c(m1_branch0,m1_branch1,m2_branch0,m2_branch1)
        #m1_b0, m1_b1, m2_b0, m2_b1
        id <- which(rb == max(rb))
        if(length(id) > 1){
          #If max is same, use min to differentiate - case arise at high levels (intermediate) with several functional cell states downstream
          minid <- which(rb == min(rb))
          if(minid == 1 | minid == 4){
            id <- id[which(!id %in% c(1,4))]
          }else{
            id <- id[which(id %in% c(1,4))]
          }
        }
        if(id == 1 | id == 4){
          #Module 1 = branch 0 & Module 2 = branch 1
          #M1 belongs to branch 0 & M2 belongs to branch 1
          #Branch 0
          dfres <- rbind(dfres,
                         data.frame("Module"="M1",
                                    "Branch"=mo[2],
                                    "Type"=typeMod[["branch0"]],
                                    "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C1),collapse = ","),
                                    "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C1)),collapse = ","),
                                    "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M1))
          #Branch 1
          dfres <- rbind(dfres,
                         data.frame("Module"="M2",
                                    "Branch"=mo[1],
                                    "Type"=typeMod[["branch1"]],
                                    "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C2),collapse = ","),
                                    "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C2)),collapse = ","),
                                    "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M2))
          
        }else if(id == 2 | id == 3){
          #Module 1 = branch 1 & Module 2 = branch 0
          #M1 belongs to branch 1 & M2 belongs to branch 0
          #Branch 0
          dfres <- rbind(dfres,
                         data.frame("Module"="M1",
                                    "Branch"=mo[1],
                                    "Type"=typeMod[["branch1"]],
                                    "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C1),collapse = ","),
                                    "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C1)),collapse = ","),
                                    "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M1))
          #Branch 1
          dfres <- rbind(dfres,
                         data.frame("Module"="M2",
                                    "Branch"=mo[2],
                                    "Type"=typeMod[["branch0"]],
                                    "TFs"=paste(names(funpartobj$cliques$Clust1[[lev[i]]]$C2),collapse = ","),
                                    "Genes"=paste(unique(unlist(funpartobj$cliques$Clust1[[lev[i]]]$C2)),collapse = ","),
                                    "Enrichment"=funpartobj$functionalenrich[[lev[i]]]$M2))
        }
      }
    }
  }else{
    print("No functional heterogeneity identified in this object")
  }
  dfres <- dfres[!is.na(dfres$Module),]
  return(dfres)
}

