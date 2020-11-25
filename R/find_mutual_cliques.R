#' Identify the best antagonistic pair of cliques
#' @param mtx single cell matrix cleaned and normalized
#' @param tf_list vector of TFs
#' @param qtarget percentile of positive targets that will appears in the gene modules
#' @param qprobInt percentiles/probs for the interactions to be considered as strong and kept for the module identification
#' @param posRatio positive ratio filter to consider a module enriched in positive 
#' @return list of pairs of negative modules
#' @author Celine Barlier
find_mutual_cliques <-
  function(mtx, tf_list, qtarget, qprobInt, posRatio){
    
#Correlation metric to infer the statistical dependencies
    m_obj = WGCNA::cor(t(mtx))
    diag(m_obj) <- 0
    
    #filter matrix to TFs only
    genes_crop = which(rownames(m_obj)%in%tf_list)
    m_obj_c = m_obj[genes_crop,genes_crop]
    
    #Keep strongest interactions
    q <- quantile(m_obj_c, probs = qprobInt)
    qp <-q[2]
    qn <- q[1]
    m_obj_c[m_obj_c < qp & m_obj_c > qn] <- 0
    
    #load matrix to graph
    g <- graph_from_adjacency_matrix(m_obj_c, mode = "undirected", weighted = T)
    #Keep only pos
    g.copy.pos = delete.edges(g,which(E(g)$weight < 0))
    
    #calculate max cliques
    clk = max_cliques(g.copy.pos, min = 3, max = 10)
    
    #If some cliques are found
    if(length(clk)!=0){
      
      #Order them
      clk = clk[order(sapply(clk,length),decreasing = T)]
      
      #unique cliques 
      co <- c()
      lc <- sapply(seq_along(clk),function(y,n,i){
        co <- c(co,n[[i]]) #already compared
        compareClk(y[[i]],clk,co)
      },y=clk, n=names(clk))
      uniques_cliques <- clk[which(lc == TRUE)]
      
      #If some unique cliques are found
      if(length(uniques_cliques)>0){
        
        #Calculate positive score for each clique - to recover macrophages results
        df_p <- foreach(x=seq(1,length(uniques_cliques)),.combine = 'c') %dopar% {
          g <- c(names(uniques_cliques[[x]]))
          m_pos <- m_obj_c[g,g]
          sum(m_pos[m_pos>0])/length(rownames(m_pos))
        }
        names(df_p) <- seq(1:length(uniques_cliques))
        df_p <- df_p[order(df_p,decreasing = T)]
        #Keep ratio positif = at least 1
        df_p <- df_p[df_p>posRatio]
        uniques_cliques <- uniques_cliques[as.numeric(names(df_p))]
        
        #If some cliques passed the first criteria
        if(length(uniques_cliques)>0){
          
          #Calculate "expression score" for each clique
          df_pexp <- foreach(y=seq(1,length(uniques_cliques)),.combine = 'c') %dopar% {
            g <- c(names(uniques_cliques[[y]]))
            m <- mtx[g,]
            mb <- m
            mb[mb>0] <- 1 #binarized
            length(which(colSums(mb) == length(g))) * (sum(m)/length(colnames(m)))
          }
          names(df_pexp) <- seq(1:length(uniques_cliques))
          #keep the clique expressed more than the mean
          df_pexp <- df_pexp[df_pexp>quantile(df_pexp,0.5)]
          uniques_cliques <- uniques_cliques[as.numeric(names(df_pexp))]
          
          #If some cliques passed the second criteria
          if(length(uniques_cliques)>1){
            
            #Calculate negative score between each pair
            df_tmp <- foreach(i=seq(1,length(uniques_cliques)),.combine='rbind',.packages = c('foreach','doParallel')) %dopar% { 
              foreach(j=seq(1,length(uniques_cliques)),.combine='rbind') %dopar% {
                #If genes in both cliques, score = 0
                if(!length(which(names(uniques_cliques[[i]]) %in% names(uniques_cliques[[j]]))) == length(names(uniques_cliques[[i]]))){
                  g <- c(names(uniques_cliques[[i]]),names(uniques_cliques[[j]]))
                  m_neg <- m_obj_c[g,g]
                  #c1exp <- which(mtx[which(rownames(mtx) %in% names(uniques_cliques[[i]])),] > 0)
                  #c1Notexp <- which(mtx[which(rownames(mtx) %in% names(uniques_cliques[[i]])),] == 0)
                  #c2exp <- which(mtx[which(rownames(mtx) %in% names(uniques_cliques[[j]])),] > 0)
                  #c2Notexp <- which(mtx[which(rownames(mtx) %in% names(uniques_cliques[[j]])),] == 0)
                  #Score
                  n <- sum(m_neg[m_neg<0]) # (sum(m_neg[m_neg<0])/length(rownames(m_neg))) (length(intersect(c1exp,c2Notexp)))/length(c1exp) * (length(intersect(c2exp,c1Notexp)))/length(c2exp)
                }else{
                  n <- 0
                }
                c(i,j,n)
              }
            }
            df_tmp <- as.data.frame(df_tmp)     
            colnames(df_tmp) <- c("C1","C2","Score")
            
            #Keep only score < 0
            df_tmp <- df_tmp[which(df_tmp$Score < 0),]
            
            #If such pattern exist continue, if not throw an empty list
            if(length(df_tmp$C1) == 0){
              return(list())
            }else{
              #Order by highest score
              df_tmp <- df_tmp[order(df_tmp$Score,decreasing = F),]
              
              #Remove all impair lines = same as pairs, A - B = B - A
              df_tmp <- df_tmp[-seq(1,length(df_tmp$Score),2),]
              
              #Frequency of each clique to be in a pair - remove the ones too frequent - most likely not unique
              #f <- table(c(df_tmp$C1,df_tmp$C2))
              #f <- sort(f,T)
              #fq1 <- as.numeric(names(f[f<quantile(f,0.25)]))
              #df_tmp <- df_tmp[which(df_tmp$C1 %in% fq1),]
              #df_tmp <- df_tmp[which(df_tmp$C2 %in% fq1),]
              
              #Keep the 10 best pairs
              l <- length(df_tmp$C1)
              len <- 10
              if(l < 10){len <- l}
              df_tmp <- df_tmp[1:len,]
              
              #List of negative pairs
              clk_neg <- list()
              for (i in seq(1,length(df_tmp$Score))) {
                clk_neg[[i]] <- list("C1"=uniques_cliques[[df_tmp$C1[i]]],"C2"=uniques_cliques[[df_tmp$C2[i]]])
              }
              
              #enrich them by targets
              #1) check target expressed in a max of same cells as clique TFs
              #2) check target almost not expressed in the antagonistic clique TFs
              #@TODO
              clk_e = sapply(clk_neg, function(x){
                list(sapply(x,function(y){
                  list(sapply(y,function(zi){
                    z = m_obj[zi,]
                    z = z[z>0]
                    z = sort(z,T)
                    z = names(z[z>quantile(z,qtarget)])
                    z = compareExp(mtx,zi,z)
                    z = names(z[z>quantile(z,qtarget)])
                    z = list(z)
                  }))
                }))
              })
              
              #mutual cliques with targets are here
              return(clk_e)
            }
            
          }else{
            return(list())
          }
          
        }else{
          return(list())
        }
        
      }else{
        return(list())
      }
      
    }else{
      return(list()) 
    }
}
