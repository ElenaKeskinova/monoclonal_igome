
# get all peptides with their frequencies
allpeps = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  V(G)$name
  
})

allfreq = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  f = V(G)$Freq
  names(f) = V(G)$name
  f
})

peps = unique(unlist(allpeps))
freqvec = lapply(peps,\(p){
  c = c(0,0,0,0)
  for(i in 1:4){
    f = allfreq[[i]][p]
    if(!is.na(f)){
      c[i] = f
    }
  }
  c
})
names(freqvec)=peps

save(freqvec,file = "allpepfreqs_perAB.RData")
# make graph
source("newgraph.R")
source("graph_to_line.R")
source("add_weights")
weight = sapply(freqvec,\(f) log2(sum(f)))
G = newgraph(peps,2)
G = set_vertex_attr(G,name = "freq_per_ab",value = frecvec)
G = set_vertex_attr(G,name = "totallogfr",value = weight)
lineg = graph_to_line(G)
w = weights(lineg,lineg$edg_lcs,v_freqs = weight)
lineg = set_edge_attr(lineg,name =  "weight",value = w)
save(lineg,file =  "big_common_graph_l.RData")

load(file =  "big_common_graph.RData")
# spectral clustering----
library(dbscan)

bc = "big_common"
lg = lineg
i = which.max(components(lg)$csize)
g = subgraph(lg,V(lg)[components(lg)$membership == i])

## Make coordinates----

arpopt=list(maxiter=100000, tol=1e-6)
L=embed_laplacian_matrix(g, no=35, which="sa", type="I-DAD", options=arpopt)#
opdim=dim_select(L$D)
coord = L$X[,2:opdim]
#proj on sphere
projS=coord/(sqrt(rowSums(coord^2)))
##dbscan----
nn = opdim*2
knnd=kNNdist(projS,nn)
knnd=sort(knnd)

result_db = chooseep_db(projS,knnd,nn)
save(result_db,file = paste0(bc,"_","dbscan.RData"))
save(coord,file = paste0(bc,"_","speccoord.RData") )
cldbsc = result_db[[1]]
ncl = result_db[[2]]
table(cldbsc)
#####################################
## make pepsets----
pepsets = gen_pepsets(lg,G,lg$edg_lcs,result_db)
freqsets = sapply(pepsets,\(set) round(weight[set]) )
logpepsets = sapply(1:length(pepsets),\(i) rep(pepsets[[i]],freqsets[[i]]))


# compare clusters to clusters of small graphs----
freqvecmat = matrix(unlist(freqvec),ncol = 4,byrow = T)
for(i in 1:4){
  print(abs[i])
  print(assortativity(G,freqvecmat[,i]))
}
ab_peps = unlist(allpepsets_w,recursive = F)
names(ab_peps) = unlist(clnames)
pep_overlap = sapply(pepsets,\(pset){
  sapply(ab_peps,\(abset){
    length(intersect(pset,abset))
  })
})
heatmap(log(pep_overlap+0.1),scale =  "none",Rowv = NA,Colv = NA)

library(universalmotif)
big_motifs = future_lapply(logpepsets,\(set)create_motif(set))
motifs_all = unlist(list(unlist(allm_w,recursive = F),bkg_m,big_motifs),recursive = F)
m_comp_b = compare_motifs(motifs_all,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
big_coord = cmdscale(m_comp_b,2)
sizes = sapply(logpepsets,length)
sizes_all = c(lengths,20000,sizes)
colors = c("red","green","blue","orange")
cols = c(colors[numcodes], "black",rep("grey",51))
plot(big_coord,cex = log10(sizes_all),col = adjustcolor(cols, alpha.f = 0.6),pch = 16,main =  "all clusters and merged graph clusters")
legend( "bottomright",legend = c(abs, "library", "merged graph"),col = c(colors, "black", "grey"),pch = 16,cex = 2)

# align peptides to make aligned ppms ----
library(parallel)
library(foreach)
require(future)
require(doParallel)
library(msa)


all_align_big = foreach(i = 3:length(logpepsets)) %dopar% {
  require(msa)
  pepset = logpepsets[[i]]
  l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  all_align_big[[i]] = l
  
}


## freq_matrices
library(future.apply)
source("freq_matrix.R")
ppms_big = future_lapply(all_align_big,\(l){
  freq_matrix(l,AA_STANDARD,ps_c = 1)
})
save(ppms_big,all_align_big,file = "ppm_big_graph.RData")


# match to epitopes ----

