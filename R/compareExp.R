#' Identify if target is expressed in similar cells as regulator and not expressed in antagonistic TFs cells
#' @param mtx single cell matrix
#' @param tf regulator name
#' @param candidatesTargets vector of potential target genes
#' @return candidateTargets vector order by the genes fullfiling all the criteria (best target -> worst target)
#' @author Celine Barlier
compareExp <- 
  function(mtx,tf,candidatesTargets){
    
    #cells for which the regulator is expressed
    reg_cells_exp <- which(mtx[tf,]>0)
    
    #cells for which the regulator is not expressed
    reg_cells_notExp <- which(mtx[tf,]<=0)
    
    tmp <- foreach(i=seq(1,length(candidatesTargets)),.combine = 'c') %dopar% {
      #gene exp cells
      g_cells_exp <- which(mtx[candidatesTargets[i],]>0)
      #gene not exp cells
      g_cells_notExp <- which(mtx[candidatesTargets[i],]==0)
      #score = product of intersection, the highest the score, the more similar the expression pattern
      (length(intersect(reg_cells_exp,g_cells_exp)) * length(intersect(reg_cells_notExp,g_cells_notExp)))
    }
    
    names(tmp) <- candidatesTargets
    tmp <- sort(tmp,T)
    
    #Order vector by score
    return(tmp)
  }