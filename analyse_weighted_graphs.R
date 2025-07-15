abs = c("Herceptin","21c","17b","b12")

library(future.apply)
library(igraph)
library(RandPro)
require(rgl)


source("add_weights.R")

for(ab in abs){
  lg = add_weights(ab)
  
  path = paste0("mixed-7graphs/",ab,"/")
  save(lg, file = paste0(path,ab,"big7-lw.RData"))
}



# spectral clustering
source("specclust.R")

for(ab in abs[1:3]){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7-lw.RData"))
  specclust(ab,lg,bc = "w")
}


# plot 3d graphics of clusters ----
ab = abs[1]
path = paste0("mixed-7graphs/",ab,"/")

load(paste0(path,ab,"big7-logw_c.RData"))

load(paste0(path,ab,"big7or.RData"))

cpl=colorRampPalette(c("#000000FF","#0050AA9F","#10AA109F","#50AF3055","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
bc = "logw_c"
load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
load(file = paste0(path,ab,"_",bc,"_","speccoord.RData"))
cldbsc = result_db[[1]]
ncl = result_db[[2]]
rndmx=form_matrix(ncol(coord),3, JLT=F)
mxsm=coord%*%rndmx
mxsm=mxsm/(sqrt(rowSums(mxsm^2)))  
plot3d(mxsm[cldbsc>0,],col=cpl(ncl-1)[cldbsc[cldbsc>0]])  
table(cldbsc)
# make pepsets 

# make pepsets ---- 
source("dbclust.R")

allpepsets_w = lapply(abs, \(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  load(paste0(path,ab,"big7-lw.RData"))
  bc = "logw"
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  pepsets = gen_pepsets(lg,G,edg_lcs,result_db)
  save(pepsets, file = paste0(path,ab,"_peps_w.RData"))
  pepsets
  
})

allfreqs = future_lapply(1:4,\(i){
  ab = abs[i]
  path = paste0("mixed-7graphs/",ab,"/")
  load(paste0(path,ab,"big7or.RData"))
  freqs= V(G)$Freq
  names(freqs) = V(G)$name
  
  sapply(allpepsets_w[[i]],\(set) freqs[set] )})

logfreqs = lapply(allfreqs,\(fr) sapply(fr,\(f) round(log2(f))))

logpepsets = lapply(1:4,\(i){
  lapply(1:length(logfreqs[[i]]),\(j){
    
    rep(allpepsets_w[[i]][[j]],logfreqs[[i]][[j]])
  })
})



clnames = sapply(1:4,\(i) sapply(1:length(logpepsets[[i]]),\(j) paste(abs[[i]],j,sep = "_")))
codes = unlist(sapply(1:4,\(i) rep(abs[i],length(logpepsets[[i]]))))
numcodes = unlist(sapply(1:4,\(i) rep(i,length(logpepsets[[i]]))))

for(i in 1:4){
  names(logpepsets[[i]]) = clnames[[i]]
}
# make motifs ----

require(universalmotif)
plan(multisession)
allm_w = future_lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(set)create_motif(set))})
allm_wp =future_lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(set)create_motif(set, pseudocount = 1))}) # with pseudocount

## aligned motifs -----

source("freq_matrix.R")
library(msa)
library(Biostrings)

all_align = lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(pepset){
    l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
    l= apply(as.matrix(l),1, paste,collapse="")
    l
    
  })
})

## freq_matrices
all_ppm = future_lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(pepset){
    require(msa)
  l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  freq_matrix(l,AA_STANDARD,ps_c = 1)
  })
  
})
save(all_ppm,all_align,file = "mixed-7graphs/ppm_pssm_w.RData")

# create pssms without and with gaps
source("freq_matrix.R")
library(msa)
library(Biostrings)
all_pssm = lapply(logpepsets,\(pepsets){
  sapply(pepsets,\(pepset){
    l=msaClustalW(AAStringSet(pepset), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
    l= apply(as.matrix(l),1, paste,collapse="")
    pssm(l)
  
  })
})

all_pssm_ = lapply(all_pssm,\(mots){
  sapply(mots,\(mot){
    mot = rbind(mot,rep(0,ncol(mot)))
    rownames(mot)[nrow(mot)] = "-"
    mot
  })
  
})

save(all_pssm,all_pssm_,all_ppm,file = "mixed-7graphs/ppm_pssm_w.RData")

# compare to clusters without weights ----
## tables with cluster distribution ----
for(ab in abs){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
  cldbsc_w = result_db[[1]]
  load(file = paste0(path,ab,"dbscan.RData"))
  cldbsc = result_db[[1]]
  print(ab)
  print(table(cldbsc,cldbsc_w))
}


# calculate kld----
kld = function(p,q){
  sum(p*log(p/q))
}
all_kld = sapply(unlist(allm_wp,recursive = F),\(mot) kld(mot@motif,bgmot))
lengths = sapply(unlist(logpepsets,recursive = F),length) 

