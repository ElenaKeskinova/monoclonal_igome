
# functions to compute paths and scores ----

longpaths = function(g){ # graph
  vs = V(g)[which(degree(g)==1)]
  peps = sapply(vs,\(v){
    paths = all_simple_paths(g,from =v,to = vs)
    peps = (sapply(paths,\(path){
      p = V(g)[names(path)]$AA
      np = length(p)
      nx = length(which(p=="-"))
      if(np-nx>4 & nx<np/3 ){ 
        p
      }
      # returns all paths as arrays
      
    }))
  })
  
}

score = function(mot,epi){ # provide a pssm and an epitope as a character array
  n = length(epi)
  if(n<=ncol(mot)){
    max(sapply(1:(ncol(mot)-n+1),\(i){
      sum(diag(mot[epi,i:(i+n-1)]))
    }))
  } 
  else{
    max(sapply(1:(n - ncol(mot)+1),\(i){
      sum(diag(mot[epi[i:(i+ncol(mot)-1)],]))
    }))
  }
}
