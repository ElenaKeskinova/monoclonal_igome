library(Biostrings)
library(igraph)
library(future.apply)
require(universalmotif)
source("newgraph.R")
source("graph_to_line.R")
source("dbclust.R")
source("rnd_graph.R")
require(parallel)
load(file = "lib_freqs.RData")
#freq = consensusMatrix(allpep)

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

load(file = "mixed-7graphs/lib_allpeps_7nC.RData")

#bkg_m = create_motif(peps_7nC)    
bgmot = freqM7
freq = rowSums(freqM7)/sum(freqM7)

n = 100 #all graphs 
ncores = detectCores()-2
sizes = sample(min(gsizes):max(gsizes),n,replace = T)

require(dbscan)
require(universalmotif)
require(igraph)
plan(multisession(workers = 20))
rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
  vc= sizes[i]
  rg = rnd_graph_reps(vc,freqs=freq,repfreqs=rep_freqs,s=i)
  
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
names(aa_prob) = names(AMINO_ACID_CODE[1:20])
aa_prob = aa_prob[order(names(aa_prob))]

n = 100 #all graphs 
ncores = detectCores()-2
sizes = sample(min(gsizes):max(gsizes),100,replace = T)

require(dbscan)
require(universalmotif)
require(igraph)
plan(multisession(workers = 20))
rnd_pepsets = future_lapply(1:n,\(i){ ### iterations for all graphs
  vc= sizes[i]
  rg = rnd_graph_reps(vc,freqs=aa_prob,repfreqs=rep_freqs,s=i+200)
  
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


rnd_align = sapply(rnd_pepsets,\(pepset){
  require(msa)
  l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(AAStringSet(l)),1, paste,collapse="")
  l
  
})

rnd_ppm = future_lapply(rnd_align,\(l){
  
  freq_matrix(l,AA_STANDARD,ps_c = 1)
})


save(rnd_align,rnd_ppm,file = "rnd_ppm_abagfreqs.RData")

