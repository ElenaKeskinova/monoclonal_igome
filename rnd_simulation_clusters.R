library(Biostrings)
library(igraph)
library(future.apply)
require(universalmotif)
source("makegraph.R")
source("graph_to_line.R")
source("dbclust.R")
source("rnd_graph.R")
require(parallel)
load(file = "lib_freqs.RData")
#freq = consensusMatrix(allpep)
freq = rowSums(freqM7)/sum(freqM7)

kld = function(p,q){
  sum(p*log(p/q))
}
abs = abs = c("Herceptin","21c","17b","b12")
gsizes = sapply(abs,\(ab){
  path = paste0("/home/elena/Documents/eli/Motifier_DataSet/mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  vcount(G)
})
load(file = "/home/elena/Documents/eli/Motifier_DataSet/mixed-7graphs/lib_allpeps_7nC.RData")

bkg_m = create_motif(peps_7nC)    
bgmot = freqM7
    
n = 100 #graphs per ab next 48 
ncores = detectCores()-2



rnd_results = lapply(gsizes[1:4],\(vc){
  
  vc = gsizes[4]
  require(dbscan)
  require(universalmotif)
  require(igraph)
  plan(multisession(workers = 16))
  res = future_lapply(1:n,\(i){ ### iterations for each antibody
    rg = rnd_graph(vc,s=i)
    
    g = rg[[1]]
    lg = rg[[2]]
    rm(rg)
    gc()
    #print(c(i,"edges" ,vcount(lg)))
    # spec clust and dbscan
    print(paste(i,"graph"))
    pepsets = dbclust_peps(g,lg)
    
    print(paste(i,"dbcl ready"))
    motifs = lapply(pepsets,\(peps) create_motif(peps,pseudocount = 1,alphabet = 'AA'))
    
    
    
    kl = sapply(motifs,\(m) kld(m@motif,bgmot))
    dist = compare_motifs(c(bkg_m,motifs),compare.to = 1,min.mean.ic = 0,max.p = 1,method = "EUC")$score
    len = sapply(pepsets,\(peps) length(peps))
    
    print(paste(i,"kl and dist"))
    return(cbind(kl,dist,len))
    
  },future.seed = T) #
  plan(sequential)
  res = sapply(1:3,\(j) unlist(sapply(1:length(res),\(i) res[[i]][,j] )))
  save(res,file = paste0("mixed-7graphs/",vc, "_res_l.RData" ))
}) #,future.seed = T

names(rnd_results) = abs[1:3]
save(rnd_results,file="mixed-7graphs/rnd_res_l.RData")

