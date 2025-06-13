
library(parallel)
library(future)
library(foreach)
library(igraph)
library(stringr)
library(doFuture)
library(msa)
library(future.apply)

source("freq_matrix.R")

abs = c("Herceptin","21c","17b","b12")

##
### make all sets of peptides for each cluster----
allpepsets = list()
plan(multisession,workers = 6)
# 
for(ab in abs){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  load(paste0(path,ab,"big7-l.RData"))
  
  i = which.max(components(lg)$csize)
  g = subgraph(lg,V(lg)[components(lg)$membership == i])
  
  load(file = paste0(path,ab,"_all_lcs.RData"))
  load(file = paste0(path,ab,"dbscan.RData"))
  cldbsc = result_db[[1]]
  ncl = result_db[[2]]
  
  ###### make pepsets
  lcssets = sapply(sort(unique(cldbsc))[2:ncl], function(i) V(g)[which(cldbsc==i)]$name)
  pepsets = future_lapply(lcssets,\(set) {
    peps = unlist(sapply(set,\(lcs) edg_lcs[which(edg_lcs[,3]==lcs),1:2]))
    unique(peps)
  })
  
  
  ic = which(igraph::components(G)$csize >=10 & igraph::components(G)$csize<max(igraph::components(G)$csize))#
  peps2 = lapply(ic, function(c) V(G)[which(igraph::components(G)$membership == c)]$name)
  
  if(length(peps2)==1){
    pepsets = append(peps2,pepsets)
  } else { pepsets = c(pepsets,peps2)}
  
  allpepsets = c(allpepsets,pepsets)
  ##
}
plan(sequential)


# create motifs ----
plan(multisession)
allp_m = future_lapply(allpepsets,\(set){ create_motif(set)})
allp_mp = sapply(pepsets, \(p){create_motif(unique(p),pseudocount = 1)}) # with pseudocount
allp_mb = append(allp_m,bkg_m) # add motif of bachground

colcodes_b = c(rep(1,39),rep(2,17),rep(3,16),rep(4,15),5)
colcodes = c(rep(1,39),rep(2,17),rep(3,16),rep(4,15))
# create pssms without and with gaps
all_pssm = lapply(allpepsets,\(pepset){
  l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  pssm(l)
})

all_pssm_ = lapply(all_pssm,\(mot){
  mot = rbind(mot,rep(0,ncol(mot)))
  rownames(mot)[nrow(mot)] = "-"
  mot
})
save(all_pssm,all_pssm_,file = "mixed-7graphs/all_pssm.RData")

## calculate distances between motifs and kld scores -----
kld = function(p,q){
  sum(p*log(p/q))
}
kl = sapply(allp_mp,\(m) kld(m@motif,bgmot)) # kld

len = sapply(allpepsets, length) # sizes of cluster
lenb = c(len, bkg_m@nsites)

# distances
m_comp_b = compare_motifs(allp_mb,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")

#pca
cmd  = cmdscale(m_comp_b,75,eig = T)
cmd_points = cmd$points[1:87,]
codes = colcodes
save(cmd_points,codes,file = "cmdresult.RData")
cluster_dist = m_comp_b
colnames(cluster_dist) = c(names(allpepsets),"bckg")
rownames(cluster_dist) = c(names(allpepsets),"bckg")
save(cluster_dist,file ="mixed-7graphs/dist_m.RData" )


