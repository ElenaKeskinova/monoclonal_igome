pepmatch = function(pep,g,simmat = simm_){ #returns all paths
  require(igraph)
  np = nchar(pep)
  pp = unlist(strsplit(pep,""))
  
  paths = lapply(1:np,\(pi){
    paa = pp[pi]
    gaa = V(g)$AA
    gaas = simmat[gaa,paa]
    good = which(gaas>mean(simmat))
    starts = V(g)[good]$name
    gaasc = gaas[good]
    if(length(starts)>0){
      if(pi == np){
        names(starts) = gaasc
        rec_paths = as.list(starts)
        
        
      } else{
        
        rec_paths = lapply(1:length(starts),\(i){
          recursive_search(starts[i],path = starts[i],pi_next = pi+1,scores = gaasc[i],pp,np,g,simmat)
        })
      }
      
      
     
    } else {
      rec_paths = as.list("gap")
      names(rec_paths) = 0
    }
     return(rec_paths)
  })
  names(paths) = 1:np
  require(purrr)
  
  paths = lapply(as.numeric(names(paths)),\(nm) {
    pnm = paths[[nm]]
    for(i in 1:(np-nm+1)){
      pnm = list_flatten(pnm)
    }
    pnm
  })
  return(paths)
}

recursive_search = function(nextv, path, pi_next,scores,pp,np,g,simmat){
  nextvs = neighbors(g, nextv)[which(!neighbors(g,nextv)$name %in% path)]
  gaas = simmat[nextvs$AA,pp[pi_next]]
  good = which(gaas>mean(simmat))
  nextvs = nextvs[good]
  gaasc = gaas[good]
  l = length(nextvs)
  pi_next = pi_next + 1
  if(l>0 & pi_next<=np){
    
    lapply(1:l,\(i){
      nextv = nextvs[i]$name
      path_ = c(path,nextv) # add chosen vertex to path
      scores_ = c(scores,gaasc[i]) # add score
      recursive_search(nextv,path_,pi_next = pi_next,scores_,pp,np,g,simmat)
    })
    
    
    
  }
  else if(l>0 & pi_next>np){ # matches exist for the final amino acid => add to existing path and return
    require(vctrs)
    for(i in 1:l){
      nextv = nextvs[i]$name
      path_ = c(path,nextv) # add chosen vertex to path
      scores_ = c(scores,gaasc[i]) # add score
      names(path_) = scores_
      return(path_)
      
    }
    
  }
  else{ # l = 0, no matches, return path up to now
    
    names(path) = scores
    return(path)
  }
}


merged_paths = function(paths,pep,distmat){ #paths is a list of lists, containing the paths with each starting letter
  require(vctrs)
  require(purrr)
  np = nchar(pep)
  pp = unlist(strsplit(pep,""))
  mp = lapply(1:(np-2),\(i){
    unique(lapply(paths[[i]],\(p){
      l = length(p)
      if(i+l+1<=np){ # path does not reach to the end of the peptide
        lapply(paths[[l+i+1]],\(p2){
          newpath = c(p,"gap",p2)
          names(newpath)[l+1] = 0
          int = intersect(p,p2)
          
          if(!("gap" %in% p) & !("gap" %in% p2)){ # both paths are not missing
            
            pathdist = distmat[p[l],p2[1]]
            if( pathdist <= 12 & length(int)==0){
              
              return(newpath)
            }
            else if(pathdist <= 12 & length(int)>=0){
              overlap = match(int,p2)
              return(newpath[1:(l+min(overlap))]) # cut up to first amino acid repetition
            }
            else{ # pathdist is > 12, take only first path
              return(p)
            }
          }
          else{
            return(p)
          }
          
        })
        
        
        
      }
      else{
        return(p)
      }
    }))
    
   
  
  })
  mp = list_drop_empty(unique(unlist(unlist(mp,recursive = F),recursive = F)))
  return(mp)
  
}

# cut paths when merging if there are overlapping letters