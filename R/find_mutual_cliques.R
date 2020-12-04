#' Identify the best antagonistic pair of cliques
#' @param mtx single cell matrix cleaned and normalized
#' @param tf_list vector of TFs
#' @param qtarget percentile of positive targets that will appears in the gene modules
#' @return list of pairs of negative modules
#' @author Celine Barlier
find_mutual_cliques <-
  function(mtx, tf_list, qtarget){
    
    #Correlation metric to infer the statistical dependencies
    m_obj = WGCNA::cor(t(mtx))
    diag(m_obj) <- 0
    
    #filter matrix to TFs only
    genes_crop = which(rownames(m_obj)%in%tf_list)
    m_obj_c = m_obj[genes_crop,genes_crop]
    
    #Keep strongest interactions
    q <- quantile(m_obj_c, probs = c(0.025,0.975))
    qp <-q[2]
    qn <- q[1]
    m_obj_c[m_obj_c < qp & m_obj_c > qn] <- 0
    
    #load matrix to graph
    g <- graph_from_adjacency_matrix(m_obj_c, mode = "undirected", weighted = T)
    #Keep only pos
    g.copy.pos = delete.edges(g,which(E(g)$weight < 0))
    
    #calculate max cliques
    clk = max_cliques(g.copy.pos, min = 3, max = 10)
    
    #If more than 2 cliques are found (needed for the antagonistic comparison)
    if(length(clk)>1){
      
      #Order them
      clk = clk[order(sapply(clk,length),decreasing = T)]
      
      #unique cliques 
      co <- c()
      lc <- foreach(i = 1:length(clk), .combine = "c", .packages = "doParallel", .export = c("compareClk")) %dopar% {
        co <- c(co,i) #already compared
        compareClk(clk[[i]],clk,co)
      }
      uniques_cliques <- clk[which(lc == TRUE)]
      
      #If some unique cliques are found
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
          #keep the top 30% cliques
          df_pexp <- df_pexp[df_pexp>quantile(df_pexp,0.7)]
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
                  #Score
                  n <- sum(m_neg[m_neg<0])/length(rownames(m_neg))
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
              
              if(length(df_tmp$C1) > 2){
                uni_c1 <- unique(df_tmp$C1)
                df_s <- foreach(i=seq(1,length(uni_c1)),.combine='rbind') %dopar% { 
                  df_unic1 <- df_tmp[which(df_tmp$C1 == uni_c1[i]),]
                  best_pair <- df_unic1[which(df_unic1$Score == min(df_unic1$Score)),]
                  c(best_pair$C1[1],best_pair$C2[1],best_pair$Score[1])
                }

                #Transformation into df
                if(length(df_s) == 1){
                  #If only 2 cliques - we got a vector
                  df_s <- data.frame("C1"=df_s[1],"C2"=df_s[2],"Score"=df_s[3])
                }else{
                  #If not we got a list of vector
                  df_s <- as.data.frame(df_s)     
                  colnames(df_s) <- c("C1","C2","Score")
                }
                
                #If few strong pairs are left (less than 10), consider all of them
                if(length(df_s[which(df_s$Score < quantile(df_s$Score,0.05)),]$C1) > 10){
                  df_s <- df_s[which(df_s$Score < quantile(df_s$Score,0.05)),]
                }
              }else{
                df_s <- df_tmp
              }

              #Consider the top 10 cliques
              df_s <- df_s[order(df_s$Score,decreasing=F),]
              l <- length(df_s$C1)
              len <- 10
              if(l < 10){len <- l}
              df_s <- df_s[1:len,]
              
              #List of negative pairs
              clk_neg <- list()
              for (i in seq(1,length(df_s$Score))) {
                clk_neg[[i]] <- list("C1"=uniques_cliques[[df_s$C1[i]]],"C2"=uniques_cliques[[df_s$C2[i]]])
              }
              
              #enrich them by targets
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
  }