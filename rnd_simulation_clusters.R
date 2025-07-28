library(Biostrings)
library(igraph)
library(future.apply)
require(universalmotif)
source("makegraph.R")
source("graph_to_line.R")
source("dbclust.R")
source("rnd_graph.R")
require(parallel)

kld = function(p,q){
  sum(p*log(p/q))
}
abs = abs = c("Herceptin","21c","17b","b12")
gsizes = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  vcount(G)
})

repfreqs = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  round(log2(V(G)$Freq))
})

rep_freqs = table(unlist(repfreqs))/sum(table(unlist(repfreqs)))



#rnd graphs for significance of clusters (kld from background)----

load(file = "lib_freqs.RData")
freq = rowSums(freqM7)/sum(freqM7)

    
n = 100 #all graphs 
ncores = detectCores()-2
sizes = sample(min(gsizes):max(gsizes),100,replace = T)

  require(dbscan)
  require(universalmotif)
  require(igraph)
  plan(multisession(workers = 20))
  rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
    vc= sizes[i]
    rg = rnd_graph(vc,freqs=freq,repfeqs=rep_freqs,s=i)
    
    g = rg[[1]]
    lg = rg[[2]]
    rm(rg)
    gc()
    #print(c(i,"edges" ,vcount(lg)))
    # spec clust and dbscan
    print(paste(i,"graph"))
    pepsets = dbclust_peps(g,lg)
    
    print(paste(i,"dbcl ready"))
    
    v_fr = V(g)$Freq
    names(v_fr) = V(g)$name
    freqs = sapply(pepsets,\(set) v_fr[set] )
    newsets = lapply(1:length(freqs),\(ii){
      rep(pepsets[[ii]],freqs[[ii]])
    })
    newsets
   
  },future.seed = T) #
  plan(sequential)
  rnd_pepsets = unlist(rnd_pepsets,recursive = F)
  save(rnd_pepsets,file = "rnd_pepsets_bkgfreqs.RData" )

  
# rnd graphs for comparison to epitopes
  
  ag_aaprob = read.csv(file = "AgAbIFprobs.csv")
  aa_prob = ag_aaprob$p
  names(aa_prob) = ag_aaprob$X
  
  n = 100 #all graphs 
  ncores = detectCores()-2
  sizes = sample(min(gsizes):max(gsizes),100,replace = T)
  
  require(dbscan)
  require(universalmotif)
  require(igraph)
  plan(multisession(workers = 20))
  rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
    vc= sizes[i]
    rg = rnd_graph(vc,freqs=aa_prob,repfeqs=rep_freqs,s=i+200)
    
    g = rg[[1]]
    lg = rg[[2]]
    rm(rg)
    gc()
    #print(c(i,"edges" ,vcount(lg)))
    # spec clust and dbscan
    print(paste(i,"graph"))
    pepsets = dbclust_peps(g,lg)
    
    print(paste(i,"dbcl ready"))
    
    v_fr = V(g)$Freq
    names(v_fr) = V(g)$name
    freqs = sapply(pepsets,\(set) v_fr[set] )
    newsets = lapply(1:length(freqs),\(ii){
      rep(pepsets[[ii]],freqs[[ii]])
    })
    newsets
    
  },future.seed = T) #
  plan(sequential)
  rnd_pepsets = unlist(rnd_pepsets,recursive = F)
  save(rnd_pepsets,file = "rnd_pepsets_abagfreqs.RData")
  
