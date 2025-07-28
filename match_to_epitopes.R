load(file = "epitope_graphs.RData")
load(file = "mixed-7graphs/ppm_pssm_logw.RData")
for(i in 1:4){
  names(all_ppm[[i]]) = clnames[[i]]
}
ppms = unlist(all_ppm,recursive = F)

require(ggseqlogo)
require(ggplot2)
require(Biostrings)


library(vctrs)
#try with blossum weights
data("BLOSUM45")
blosum = BLOSUM45[1:20,]
blosum = cbind(blosum,rep(0,20))
colnames(blosum)[ncol(blosum)] = "-"

max_scores_all = future_sapply(1:4,\(ab_i){
  paths = longpaths(g_full[[ab_i]])
  sapply(ppms,\(mot){
    
    scores = unlist(sapply(paths,\(paths2) {
      paths2 = list_drop_empty(paths2)
      sapply(paths2,\(path){
        n = length(path)
        if(n<=ncol(mot)){
          max(sapply(1:(ncol(mot)-n+1),\(i){
            #sum(diag(mot[epi,i:(i+n-1)]))
            m = mot[,i:(i+n-1)]*cbind(blosum[,path]) # part of motif times substitution frequency of path aas
            mean(colSums(m))
          }))
        } 
        else{
          max(sapply(1:(n - ncol(mot)+1),\(i){
            #sum(diag(mot[epi[i:(i+ncol(mot)-1)],]))
            m = mot*cbind(blosum[,path[i:(i+ncol(mot)-1)]])
            mean(colSums(m))
          }))
        }
      })
    }))
    max(scores)
  })
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

# all scores
allpaths = lapply(g_full,\(g){
  paths = longpaths(g)
  paths = unlist(paths,recursive = F)
  paths = list_drop_empty(paths)
  sapply(paths,\(p) paste0(p,collapse = ""))
})

scores_all = future_sapply(1:4,\(ab_i){
  paths = longpaths(g_full[[ab_i]])
  sapply(ppms,\(mot){
    paths = unlist(paths,recursive = F)
    paths = list_drop_empty(paths)
    scores = unlist(sapply(paths,\(path) {
      
        n = length(path)
        if(n<=ncol(mot)){
          max(sapply(1:(ncol(mot)-n+1),\(i){
            #sum(diag(mot[epi,i:(i+n-1)]))
            
            m = mot[,i:(i+n-1)]*cbind(blosum[,path]) # part of motif times substitution frequency of path aas
            mean(colSums(m))
          }))
        } 
        else{
          max(sapply(1:(n - ncol(mot)+1),\(i){
            #sum(diag(mot[epi[i:(i+ncol(mot)-1)],]))
            m = mot*cbind(blosum[,path[i:(i+ncol(mot)-1)]])
            mean(colSums(m))
          }))
        }
    
    }))
    scores
  })
})

for(i in 1:4){
  rownames(scores_all[[i]]) = allpaths[[i]]
}

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


# random motifs ----

ag_aaprob = read.csv(file = "AgAbIFprobs.csv")
aa_prob = ag_aaprob$p
names(aa_prob) = ag_aaprob$X
N = sample(20:1000,200,replace = T)
rnd_mot = sapply(N,\(n){
  peps = sapply(1:n,\(i){
    paste0(sample(AA_STANDARD,7,prob = aa_prob,replace = T),collapse = "")
  })
  l=msaClustalW(AAStringSet(peps), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  freq_matrix(l,AA_STANDARD,ps_c = 1)

})

rnd_scores = future_sapply(allpaths,\(paths){
  sapply(rnd_mot,\(mot){
    scores = unlist(sapply(paths,\(path) {
      path = unlist(strsplit(path,""))
      score_blosum(mot,path,blosum)
      
    }))
    
    max(scores)
    
  })
})
save(rnd_scores, file = "rnd_motif_scores.RData")
hist(rnd_scores)


##rnd motifs which are from similar peptides----
N = sample(20:1000,200,replace = T)
rnd_good = sapply(N,\(n){
  peps = list()
  peps[[1]] = sample(AA_STANDARD,7,prob = aa_prob,replace = T)
  for(i in 2:n){
    j = sample(1:(i-1),1)
    k = sample(1:7,2)
    aa = sample(AA_STANDARD,2,prob = aa_prob,replace = T)
    pep = peps[[j]]
    pep[k] = aa 
    peps[[i]] = pep
  }
  peps = sapply(peps,\(pep) paste0(pep,collapse = ""))
  l=msaClustalW(AAStringSet(peps), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  freq_matrix(l,AA_STANDARD,ps_c = 1)
})

save(rnd_good,rnd_mot,file = "rnd_motifs.RData")

rnd_scores_good = future_sapply(allpaths,\(paths){
  sapply(rnd_good,\(mot){
    scores = unlist(sapply(paths,\(path) {
      path = unlist(strsplit(path,""))
      score_blosum(mot,path,blosum)
      
    }))
    
    max(scores)
    
  })
})
save(rnd_scores,rnd_scores_good, file = "rnd_motif_scores.RData")
hist(rnd_scores_good)

## also with reps----
repfreqs = table(unlist(logfreqs))/sum(table(unlist(logfreqs))) 
N = sample(10:1000,200,replace = T)
plan(multisession)
rnd_mot_reps = future_sapply(N,\(n){
  reps = sample(names(repfreqs),n,prob = repfreqs,replace=T)
  peps = list()
  peps[[1]] = sample(AA_STANDARD,7,prob = aa_prob,replace = T)
  for(i in 2:n){
    j = sample(1:(i-1),1)
    k = sample(1:7,2)
    aa = sample(AA_STANDARD,2,prob = aa_prob,replace = T)
    pep = peps[[j]]
    pep[k] = aa 
    peps[[i]] = pep
  }
  peps = sapply(peps,\(pep) paste0(pep,collapse = ""))
  peps = rep(peps,reps)
  l=msaClustalW(AAStringSet(peps), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  freq_matrix(l,AA_STANDARD,ps_c = 1)
},future.seed = T)

save(rnd_good,rnd_mot,rnd_mot_reps,file = "rnd_motifs.RData")

rnd_scores_reps = future_sapply(allpaths,\(paths){
  sapply(rnd_mot_reps,\(mot){
    scores = unlist(sapply(paths,\(path) {
      path = unlist(strsplit(path,""))
      score_blosum(mot,path,blosum)
      
    }))
    
    max(scores)
    
  })
})
save(rnd_scores,rnd_scores_good,rnd_scores_reps, file = "rnd_motif_scores.RData")
hist(rnd_scores_reps)


# load from big computer
load(file ='rnd_motifs_.RData')
load(file = "rnd_motif_scores_.RData")

#plot with random scores as controls----
for(i in 1:4){
  l = sapply(logfreqs[[i]],sum)
  o = order(max_scores_all[which(codes==abs[i]),i])
  sc = max_scores_all[which(codes==abs[i]),i][o]
  control =rnd_scores_reps[,i]
  e = ecdf(control)
  qt = e(sc)
  p = 1-qt
  p = p.adjust(p)
  
  qth = which(p< 0.1)
  
  names(p) = allnames[which(codes==abs[i])][o]
  print(p)
  
  plot(sc, col = colors[i], ylim = range(max_scores_all),pch = 16,cex = 1.3,main = paste(abs[i], "with rnd motifs"),ylab = "scores",xlab = "smallest to largest cluster")
  
  colin = "lightgrey"
  colbd = "grey"
  vioplot(control,add = T,at = 1.5,wex = 1.5, col=colin, border=colbd, rectCol=colbd, lineCol=colbd)
  points(sc, col = colors[i],pch = 16,cex = 1.3)
  points(sc, cex = 1.3)
  points(1,max_epi_scores[i,i],cex = 1.3)
  if(min(p)<0.1){
    x = (1:length(sc)+0.2)
    y = (sc+ sd(max_scores_all)/4)
    text(x[qth], y[qth], labels = round(p[qth],2), cex = 1.2)
    
  }
  
}
allal = unlist(all_align,recursive = F)

pdf(file = "all_mot_dbscan.pdf",width=5, height = 5)

for (n in allnames[1:93]){
  logoal(allal[[n]],nm = paste(n,"peps:",length(allal[[n]])))
}
dev.off()

epi_p = sapply(1:4,\(i){
  l = sapply(logfreqs[[i]],sum)
  o = order(max_scores_all[which(codes==abs[i]),i])
  sc = max_scores_all[which(codes==abs[i]),i][o]
  control =rnd_scores_reps[,i]
  e = ecdf(control)
  qt = e(sc)
  p = 1-qt
  p = p.adjust(p)
  
  qth = which(p< 0.1)
  
  names(p) = allnames[which(codes==abs[i])][o]
  p
})
