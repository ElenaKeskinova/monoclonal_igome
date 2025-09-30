load(file = "epitope_graphs.RData")
load(file = "mixed-7graphs/ppm_pssm_logw.RData")
source("paths_scores_functions.R")

clnames = sapply(1:4,\(i) sapply(1:length(logpepsets[[i]]),\(j) paste(abs[[i]],j,sep = "_")))
ppms = unlist(all_ppm,recursive = F)
names(ppms) = unlist(clnames)

ppms = lapply(ppms,\(m){
  nc = ncol(m)
  m = rbind(m,rep(0,nc))
  rownames(m)[21] = "-"
  m
})

pssms = unlist(pssm_7_clean2,recursive = F)

require(ggseqlogo)
require(ggplot2)
require(Biostrings)



#try with blossum weights
data("BLOSUM45")
blosum = BLOSUM45[1:20,]
blosum = cbind(blosum,rep(0,20))
colnames(blosum)[ncol(blosum)] = "-"




# all scores
library(vctrs)
allpaths = lapply(g_full,\(g){
  paths = longpaths(g)
  paths = unlist(paths,recursive = F)
  paths
})

scores_all = future_sapply(1:4,\(ab_i){
  paths = allpaths[[ab_i]]
  sapply(pssms,\(mot){
    
    scores = unlist(sapply(paths,\(path) {
        path = unlist(strsplit(path,""))
        score(mot,path)
    }))
    scores
  })
})

for(i in 1:4){
  rownames(scores_all[[i]]) = allpaths[[i]]
}

max_scores_all = future_sapply(scores_all,\(mat){
  apply(mat,2,max)
})

hist(max_scores_all)
# motifs from each graph vs all paths -> access best scores----
gmots = lapply(allpaths,\(paths){
  l=msaClustalW(AAStringSet(paths), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  freq_matrix(l,AA_STANDARD,ps_c = 1)

})
names(gmots) = abs
##scores ----
epi_scores = future_lapply(allpaths,\(paths){
  
  sapply(gmots,\(mot){
    
      scores = unlist(sapply(paths,\(path) {
        path = unlist(strsplit(path,""))
        score_blosum(mot,path,blosum)
        
      }))
      names(scores) = paths
      scores

  })
})

max_epi_scores = future_sapply(epi_scores,\(mat){
  apply(mat,2,max)
})


# plot score distributions
library(vioplot)
colors = c("red","green","blue","orange")
for(i in 1:4){
  l = sapply(logfreqs[[i]],sum)
  o = order(max_scores_all[which(codes==abs[i]),i])
  sc = max_scores_all[which(codes==abs[i]),i][o]
  control = max_scores_all[,i]
  e = ecdf(control)
  qt = e(sc)
  qth = which(qt>0.9)
  plot(sc, col = colors[i], ylim = range(max_scores_all),pch = 16,cex = 1.3,main = abs[i],ylab = "scores",xlab = "smallest to largest cluster")
  
  colin = "lightgrey"
  colbd = "grey"
  vioplot(control,add = T,at = 1.5,wex = 1.5, col=colin, border=colbd, rectCol=colbd, lineCol=colbd)
  points(sc, col = colors[i],pch = 16,cex = 1.3)
  points(sc, cex = 1.3)
  if(max(qt)>0.9){
    x = (1:length(sc)+0.2)
    y = (sc+ sd(max_scores_all)/4)
    text(x[qth], y[qth], labels = round(qt[qth],2), cex = 1.2)
    
  }
  
}
# rnd scores on big computer ----
# scores with pssm from background frequencies ----
plan(multisession(workers = 20))
rnd_scores_bkg = future_sapply(allpaths,\(paths){
  sapply(rnd_pssm_bkg,\(mot){
    scores = unlist(sapply(paths,\(path) {
      path = unlist(strsplit(path,""))
      score(mot,path)
      
    }))
    
    max(scores)
    
  })
})
save(rnd_scores_bkg, file = "rnd_motif_scores_bkg.RData")

# load from big computer
load(file ='rnd_motifs_.RData')
load(file = "rnd_motif_scores_.RData")
load(file = "rnd_motif_scores_bkg_2.RData") # with bkg frequencies, from pssm

hist(sample(rnd_scores_bkg,500))
hist(max_scores_all, col = adjustcolor("green", alpha = 0.5),add = T)

 #plot with random scores as controls----

for(i in 1:4){
  l = sapply(logfreqs[[i]],sum)
  o = order(max_scores_all[which(codes==abs[i]),i])
  sc = max_scores_all[which(codes==abs[i]),i][o]
  control =rnd_scores_bkg[,i]
  e = ecdf(control)
  qt = e(sc)
  p = 1-qt
  p = p.adjust(p)
  
  qth = which(p< 0.1)
  
  names(p) = clnames[[i]][o]
  print(p)
  
  plot(sc, col = colors[i], ylim = range(max_scores_all),pch = 16,cex = 1.3,main = paste(abs[i], "with rnd motifs"),ylab = "scores")
  
  colin = "lightgrey"
  colbd = "grey"
  vioplot(control,add = T,at = 1.5,wex = 1.5, col=colin, border=colbd, rectCol=colbd, lineCol=colbd)
  points(sc, col = colors[i],pch = 16,cex = 1.3)
  points(sc, cex = 1.3)
  # points(1,max_epi_scores[i,i],cex = 1.3)
  if(min(p)<0.1){
    x = (1:length(sc)+0.2)
    y = (sc+ sd(max_scores_all)/4)
    text(x[qth], y[qth], labels = round(p[qth],2), cex = 1.2)
    
  }
  
}

# plot on one plot

plot(1, type = "n", xlab = "", ylab = "similarity score",ylim = range(max_scores_all),xlim = c(0,8.5),axes = F, frame.plot = FALSE)
#axis(1, labels = FALSE)  # x-axis without labels
axis(2, labels = T)  
for(i in 1:4){
  
  sc = max_scores_all[which(codes==abs[i]),i]
  control =rnd_scores_bkg[,i]
  e = ecdf(control)
  qt = e(sc)
  p = 1-qt
  p = p.adjust(p)
  
  qth = which(p< 0.05)
  
  names(p) = clnames[[i]]
  
  pos = i*2
  points(rep(pos,length(sc)),sc, col = colors[i],pch = 16,cex = 1.3,main = paste(abs[i], "with rnd motifs"),ylab = "scores")
  
  colin = "lightgrey"
  colbd = "grey"
  at = i*2-1
  vioplot(control,add = T,at = at,wex = 1.5, col=colin, border=colbd, rectCol=colbd, lineCol=colbd)
  
  # points(1,max_epi_scores[i,i],cex = 1.3)
  if(min(p)<0.05){
    x = (i*2+0.2)
    y = (sc+ sd(max_scores_all)/4)
    text(x, y[qth], labels = round(p[qth],2), cex = 1.2)
    
  }
  
}


allal = unlist(all_align,recursive = F)

pdf(file = "all_mot_dbscan.pdf",width=5, height = 5)

for (n in allnames[1:93]){
  logoal(allal[[n]],nm = paste(n,"peps:",length(allal[[n]])))
}
dev.off()

# calculate p values ----
epi_p = sapply(1:4,\(i){
  l = sapply(logfreqs[[i]],sum)
  o = order(max_scores_all[which(codes==abs[i]),i])
  sc = max_scores_all[which(codes==abs[i]),i][o]
  control =rnd_scores_bkg[,i]
  e = ecdf(control)
  qt = e(sc)
  p = 1-qt
  p = p.adjust(p)
  
  qth = which(p< 0.1)
  
  names(p) = allnames[which(codes==abs[i])][o]
  p
})

# choose peptides from good clusters to test with docking ----
good_cl = sapply(epi_p,\(ps){
  best = names(sort(ps))[1:2]
})
## load graphs
graphs = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"big7or.RData"))
  G
})
## subgraphs from clusters
subg_all = lapply(1:4,\(abi){
  sapply(allpepsets_wm[[abi]],\(set){
    
    subgraph(graphs[[abi]], set)
  })
})
for(i in 1:4){ names(subg_all[[i]]) = names(allpepsets_wm[[i]])}
subg_good = lapply(1:4,\(abi){
  sapply(good_cl[,abi],\(name){
    subg_all[[abi]][[name]]
  })
})
names(subg_good) = abs
sapply(unlist(subg_good, recursive = F), plot.igraph)
sapply(unlist(subg_good, recursive = F),\(g) components(g)$csize)

subg_centrality = lapply(subg_all,\(gab){
  lapply(gab,\(g) eigen_centrality(g)$vector)
})
subg_central = lapply(subg_centrality,\(gab){
  lapply(gab,\(v) names(which.max(v)))
})

bestmatch = sapply((epi_p),\(ps) names(which.min(ps)))
worstmatch = sapply((epi_p),\(ps) names(which.max(ps)))

clscores = lapply(1:4,\(i){
  sc = max_scores_all[which(codes==abs[i]),i]
  control =rnd_scores_bkg[,i]
  e = ecdf(control)
  qt = e(sc)
  names(qt) = names(sc)
  qt
})

scoresmimo_to_clust = lapply(1:4,\(i){
  lapply(1:length(allpepsets_w[[i]]),\(j){
    mot = pssm_7_clean[[i]][[j]]
    set = allpepsets_w[[i]][[j]]
    sapply(set,\(pep){
      
        score(mot,unlist(strsplit(pep,"")))
      
    })
  })
})

for(i in 1:4){
  for(j in 1:length(allpepsets_w[[i]])){
    set = allpepsets_w[[i]][[j]]
    plot(subg_centrality[[i]][[j]][set],scoresmimo_to_clust[[i]][[j]][set], xlab = "centrality", ylab = "score", main = paste(ab[i],"cluster",j))
  }
}

# try scores with number of aa with pssm>1
plan(multisession(workers = 4))
scores_all2 = future_sapply(1:4,\(ab_i){
  paths = allpaths[[ab_i]]
  sapply(pssms,\(mot){
    
    scores = unlist(sapply(paths,\(path) {
      path = unlist(strsplit(path,""))
      score_binary(mot,path)
    }))
    names(scores) = paths
    scores
  })
})

for(i in 1:4){
  rownames(scores_all2[[i]]) = allpaths[[i]]
}

max_scores_all2 = future_sapply(scores_all2,\(mat){
  apply(mat,2,max)
})

hist(max_scores_all2)

scores_per_ab = sapply(1:4,\(i){
  o = order(max_scores_all[which(codes==abs[i]),i])
  sc = max_scores_all[which(codes==abs[i]),i][o]
  sc
})
scores_per_ab2 = sapply(1:4,\(i){
  o = order(max_scores_all2[which(codes==abs[i]),i])
  sc = max_scores_all2[which(codes==abs[i]),i][o]
  sc
})

for(i in 1:4){
  c = cor(scores_per_ab[[i]],scores_per_ab2[[i]][names(scores_per_ab[[i]])])
  print(c)
}
for(i in 1:4){
  plot(scores_per_ab[[i]],scores_per_ab2[[i]][names(scores_per_ab[[i]])], main = abs[[i]])
  
}
