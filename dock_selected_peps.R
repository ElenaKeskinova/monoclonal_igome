# best peps from matching to epitopes and closest to the background

peps_best = unlist(subg_central)[bestmatch]
peps_closest = unlist(subg_central)[closest]
peps_worst = unlist(subg_central)[worstmatch]

sapply(c(peps_best, peps_closest, peps_worst), print)
# get epitope contacts----
# all contacts between epitope and paratope
abs = c("Herceptin", "21c"   ,  "17b",       "b12" )

epi_contacts = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/",ab,"_contacts_epi.txt"))
})
epi_contacts = lapply(epi_contacts,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ",")); sapply(aas, \(aa) substr(aa,1,nchar(aa)-1))})
  ab
})
names(epi_contacts) = abs


# get model contacts ----

allfiles = list.files("dock_motifier_mimotopes/best_haddock_models", pattern = ".\\.txt")

models_best_contacts = lapply(abs, \(ab){
  lapply(1:2,\(i){ # 2 top clusters
    f = allfiles[grep(paste0(ab,"_best_top",i),allfiles)]
    
    lapply(f, \(f1){# 4 models per cluster
      file = paste0("dock_motifier_mimotopes/best_haddock_models/",f1)
      read_cont(file) 
    })
    
  })
})
names(models_best_contacts) = abs

models_closest_contacts = lapply(abs, \(ab){
  lapply(1:2,\(i){ # 2 top clusters
    f = allfiles[grep(paste0(ab,"_closest_top",i),allfiles)]
    
    lapply(f, \(f1){# 4 models per cluster
      file = paste0("dock_motifier_mimotopes/best_haddock_models/",f1)
      read_cont(file) 
    })
    
  })
})
names(models_closest_contacts) = abs

models_worst_contacts = lapply(abs[1:2], \(ab){
  lapply(1:2,\(i){ # 2 top clusters
    f = allfiles[grep(paste0(ab,"_worst_top",i),allfiles)]
    
    lapply(f, \(f1){# 4 models per cluster
      file = paste0("dock_motifier_mimotopes/best_haddock_models/",f1)
      read_cont(file) 
    })
    
  })
})
names(models_worst_contacts) = abs[1:2]

# contact matrix ----

cont_mats_best = lapply(1:4,\(i){
  mimo = peps_best[i]
  epi_c = epi_contacts[[i]]
  lapply(1:2,\(j){ # 2 clusters
    cont_mat(epi_c, models_best_contacts[[i]][[j]], mimo)
  })
})
cont_mats_closest = lapply(1:4,\(i){
  mimo = peps_closest[i]
  epi_c = epi_contacts[[i]]
  lapply(1:2,\(j){ # 2 clusters
    cont_mat(epi_c, models_closest_contacts[[i]][[j]], mimo)
  })
})
names(cont_mats_best) = abs
names(cont_mats_closest) = abs

cont_mats_worst = lapply(1:2,\(i){
  mimo = peps_worst[i]
  epi_c = epi_contacts[[i]]
  lapply(1:2,\(j){ # 2 clusters
    cont_mat(epi_c, models_worst_contacts[[i]][[j]], mimo)
  })
})
names(cont_mats_worst) = abs[1:2]

lapply(cont_mats_best,\(mat) lapply(mat, colSums))
lapply(cont_mats_closest,\(mat) lapply(mat, colSums))

# make metrics: % contacts of mimo in epi and of epi in mimo
overlap_percent_best = lapply(1:4,\(i){
  
  epi_c = epi_contacts[[i]]
  epi_c_all = unique(unlist(epi_c$AB_contacts))
  nepic = length(epi_c_all)
  sapply(1:2,\(j){ # 2 clusters
    cl4models = models_best_contacts[[i]][[j]]
    mimo_c_all = lapply(cl4models,\(df) unique(unlist(df$AB_contacts)))
    nmimoc = sapply(mimo_c_all, length)
    ofmimo = sapply(1:4,\(k) length(intersect(epi_c_all, mimo_c_all[[k]]))/nmimoc[k])
    ofepi = sapply(1:4,\(k) length(intersect(epi_c_all, mimo_c_all[[k]]))/nepic)
    result = c(mean(ofmimo),mean(ofepi))
    names(result) = c("ofmimo","ofepi")
    return(result)
  })
})


overlap_percent_closest = lapply(1:4,\(i){
  
  epi_c = epi_contacts[[i]]
  epi_c_all = unique(unlist(epi_c$AB_contacts))
  nepic = length(epi_c_all)
  sapply(1:2,\(j){ # 2 clusters
    cl4models = models_closest_contacts[[i]][[j]]
    mimo_c_all = lapply(cl4models,\(df) unique(unlist(df$AB_contacts)))
    nmimoc = sapply(mimo_c_all, length)
    ofmimo = sapply(1:4,\(k) length(intersect(epi_c_all, mimo_c_all[[k]]))/nmimoc[k])
    ofepi = sapply(1:4,\(k) length(intersect(epi_c_all, mimo_c_all[[k]]))/nepic)
    result = c(mean(ofmimo),mean(ofepi))
    names(result) = c("ofmimo","ofepi")
    return(result)
  })
})

overlap_percent_worst = lapply(1:2,\(i){
  
  epi_c = epi_contacts[[i]]
  epi_c_all = unique(unlist(epi_c$AB_contacts))
  nepic = length(epi_c_all)
  sapply(1:2,\(j){ # 2 clusters
    cl4models = models_worst_contacts[[i]][[j]]
    mimo_c_all = lapply(cl4models,\(df) unique(unlist(df$AB_contacts)))
    nmimoc = sapply(mimo_c_all, length)
    ofmimo = sapply(1:4,\(k) length(intersect(epi_c_all, mimo_c_all[[k]]))/nmimoc[k])
    ofepi = sapply(1:4,\(k) length(intersect(epi_c_all, mimo_c_all[[k]]))/nepic)
    result = c(mean(ofmimo),mean(ofepi))
    names(result) = c("ofmimo","ofepi")
    return(result)
  })
})

plot(c(NA,NA),xlim = c(0,1),ylim = c(0,1), xlab = "% of mimo", ylab = "% of epi")
for(res in overlap_percent_best){
  points(x = res["ofmimo",],y = res["ofepi",],pch = 16, col = c("blue","green"))
}
for(res in overlap_percent_closest){
  points(x = res["ofmimo",],y = res["ofepi",],pch = 16, col = c("red","orange"))
}


for(i in 1:4){
  ab  = abs[i]
  resb = overlap_percent_best[[i]]
  resc = overlap_percent_closest[[i]]
  plot(x = resb["ofmimo",],y = resb["ofepi",],pch = 16, col = c("blue","green"), main = ab, xlim = c(0,1),ylim = c(0,1))
  points(x = resc["ofmimo",],y = resc["ofepi",],pch = 16, col = "red")
  if(i<3){
    recw = overlap_percent_worst[[i]]
    points(x = recw["ofmimo",],y = recw["ofepi",],pch = 16, col = "orange")
  }
}

# plot of mimo vs scores

scoresbest = unlist(clscores)[bestmatch]
scoresclose = unlist(clscores)[closest]
scoresworst = unlist(clscores)[worstmatch]

for(i in 1:4){
  ab  = abs[i]
  resb = overlap_percent_best[[i]]
  resc = overlap_percent_closest[[i]]
  plot(x = resb["ofmimo",],y = rep(scoresbest[i],2),pch = 16,cex = 2, col = "green", main = ab, xlim = c(0,1),ylim = c(0,1))
  points(x = resc["ofmimo",],y = rep(scoresclose[i],2),pch = 16,cex = 2, col = "red")
  if(i<3){
    recw = overlap_percent_worst[[i]]
    points(x = recw["ofmimo",],y = rep(scoresworst[i],2),pch = 16,cex = 2, col = "orange")
  }
}


# plot epi graphs with size being number of overlapping contacts
twodim_scale = lapply(dist_mx,\(dmat){
  
  s = cmdscale(dmat,k=2)
  rownames(s) = rownames(dmat)
  s
})

g_coords = lapply(1:4,\(i){
  gf = g_full[[i]]
  coords1 = twodim_scale[[i]]
  vx = V(gf)$name[which(V(gf)$AA=="-")]
  coords2 = sapply(vx,\(x){
    nb = neighbors(gf,x)$name
    nbc = coords1[nb,]
    xc = c(mean(nbc[,1]),mean(nbc[,2]))
    xc
  } )
  rbind(coords1,t(coords2))
})

scales = lapply(cont_mats_best,\(model){
  lapply(model,\(clusterm){
    intensity = rowSums(clusterm)
    position = sapply(1:nrow(clusterm),\(j) sum(clusterm[j,]*(1:ncol(clusterm)))/sum(clusterm[j,]))
    mat = cbind(intensity,position)
    rownames(mat) = rownames(clusterm)
    mat
  })
})

# save picture of graph with size indicating number of contacts
for(i in 1:4){
  gf = g_full[[i]]
  scl = scales[[i]][[1]][,1] # overall common contacts
  epicnt = sapply(epi_contacts[[i]]$AB_contacts,length)
  names(epicnt) = epi_contacts[[i]]$AG_contact
  # number of contacts of each aa from epitope, for normalisation of common contacts
  epicnt = epicnt[names(scl)]
  #pslog = pseudo_log_trans(sigma = 4)
  #coord = pslog$transform(g_coords[[i]][V(gf)$name,])
  coord = (g_coords[[i]][V(gf)$name,])
  sizes = rep(1, vcount(gf))
  
  names(sizes) = V(gf)$name
  sizes[names(scl)] =(scl + 1)/epicnt
  png(paste0("dock_motifier_mimotopes/",abs[[i]],"_epi_graph_bestpep.png"),width = 800,height = 800)
  plot(gf, vertex.size = sizes,main = "", vertex.label =V(gf)$AA,layout = coord,vertex.frame.color = NA, vertex.label.cex = 3.2, main.size = 2)
  text(0.8,1, paste("Mimotope:",peps_best[i],"\nBest path:",bestpaths[i]), cex = 1.8 )
  dev.off()
  }

for(i in 1:1){
  gf = g_full[[i]]
  scl = scales[[i]][[1]][,1] # overall common contacts
  epicnt = sapply(epi_contacts[[i]]$AB_contacts,length) # contacts per aa of epitope
  names(epicnt) = epi_contacts[[i]]$AG_contact
  # number of contacts of each aa from epitope, for normalisation of common contacts
  epicnt = epicnt[names(scl)]
  #pslog = pseudo_log_trans(sigma = 4)
  #coord = pslog$transform(g_coords[[i]][V(gf)$name,])
  coord = (g_coords[[i]][V(gf)$name,])
  sizes = rep(0.5, vcount(gf))
  
  names(sizes) = V(gf)$name
  sizes[names(scl)] =(scl + 1)/epicnt
  png(paste0("dock_motifier_mimotopes/",abs[[i]],"_epi_graph_bestpep_legend.png"),width = 800,height = 800)
  plot(gf, vertex.size = sizes,main = "", vertex.label =V(gf)$AA,layout = coord,vertex.frame.color = NA, vertex.label.cex = 3.2, main.size = 2)
  text((max(coord[,1])-10),(max(coord[,2])-10), paste("Mimotope:",peps_best[i],"\nBest path:",bestpaths[i]), cex = 3 )
  text(1,1, paste("Mimotope:",peps_best[i],"\nBest path:",bestpaths[i]), cex = 3 )
  
  sizesl <- c(0.5, 1, 2)  # example sizes to display
  legend("bottomright", legend = paste(c(50,100,100),"% contact overlap with",c(1,1,2),"\n AA of mimotope"),
       pt.cex = sizesl*4,  # scale to match node size
       pch = 16, col = "orange",
       title = "Node size")

  dev.off()
}

# no sizes
for(i in 1:1){
  gf = g_full[[i]]
  scl = scales[[i]][[1]][,1] # overall common contacts
  epicnt = sapply(epi_contacts[[i]]$AB_contacts,length) # contacts per aa of epitope
  names(epicnt) = epi_contacts[[i]]$AG_contact
  # number of contacts of each aa from epitope, for normalisation of common contacts
  epicnt = epicnt[names(scl)]
  #pslog = pseudo_log_trans(sigma = 4)
  #coord = pslog$transform(g_coords[[i]][V(gf)$name,])
  coord = (g_coords[[i]][V(gf)$name,])
  sizes = rep(0.5, vcount(gf))
  
  names(sizes) = V(gf)$name
  sizes[names(scl)] =(scl + 1)/epicnt
  
  plot(gf, main = "", vertex.label =V(gf)$AA,layout = coord,vertex.frame.color = NA, vertex.label.cex = 2, main.size = 2)
  #text((max(coord[,1])-10),(max(coord[,2])-10), paste("Mimotope:",peps_best[i],"\nBest path:",bestpaths[i]), cex = 3 )
  #text(1,1, paste("Mimotope:",peps_best[i],"\nBest path:",bestpaths[i]), cex = 3 )
  
  
}


bestpaths = sapply(1:4,\(i){
  nm = names(peps_best)[i]
  path = which.max(scores_all[[i]][,nm])
  epi = unlist(strsplit(names(path),""))
  n = length(epi)
  mot = pssm_7_clean[[i]][[nm]]
  
  if(n<=ncol(mot)){
    sc=sapply(1:(ncol(mot)-n+1),\(i){
      d = (diag(mot[epi,i:(i+n-1)]))
      sum(d[which(d>0)])
    })
  } else{
    sc=sapply(1:(n - ncol(mot)+1),\(i){
      d = (diag(mot[epi[i:(i+ncol(mot)-1)],]))
      sum(d[which(d>0)])
    })
  }
  
  j = which.max(sc)
  substr(names(path),j,(j+6))
  
})
