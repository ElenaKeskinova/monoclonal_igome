library(Biostrings)
source("paths_scores_functions.R")
# use all_pssm_ from the script create_motifs.R
# take files with results from Pymol for contact amino acids on antigen surfaces, 
# obtained from structures 1yyl,3lqa,2ny7,1n8z from pdb structure with contact cutoff 5A
# make graphs from epitope amino acids, introduce gap nodes where distance is over 6 A , take all paths and find highest overlap score to each cluster


AA_code = sapply(AMINO_ACID_CODE,toupper)
AA_code = c(AA_code,"XXX")
names(AA_code)[27] = "-"
abs = c("Herceptin","21c","17b","b12")

distances = lapply(abs,\(ab){ # distances between amino acids closer than 5 A to antibody, from PDB, calculated with Pymol
  df = read.table(paste0("mixed-7graphs/epi_",ab,".txt"),header = T)
  df$aan1 = sapply(df$Residue1, substr, 1,3)
  df$aan2 = sapply(df$Residue2, substr, 1,3)
  df[which(df$aan1 %in% AA_code & df$aan2 %in% AA_code),] # only valid amino acids
  
})

# create distance matrix

dist_mx = lapply(distances,\(df){
  alaa = unique(unlist(df[,1:2]))
  m = matrix(0, nrow=length(alaa), ncol=length(alaa),
             dimnames=list(alaa, alaa))
  for(i in 1:nrow(df)){
    e1 = df[i,1]
    e2 = df[i,2]
    d = df[i,3]
    m[e1,e2] = d
    m[e2,e1] = d
  }
  m
})

# make graphs ----
### minimal graphs ----
# start with small graphs and the add edges to connect all amino acids of the epitopes

d = 6
ds = lapply(distances,\(df){
  df$Distance<=d
  
})

graphs_min = lapply(1:4,\(i){
  require(igraph)
  e = ds[[i]]
  edges = cbind(distances[[i]]$Residue1,distances[[i]]$Residue2)[e,]
  g = graph_from_edgelist(edges, directed = F)
  aas = sapply(V(g)$name,\(nm){
    c3 = substr(nm,1,3) 
    names(AA_code[which(AA_code==c3)])
    
  } )
  set_vertex_attr(g,name = "AA", value = unlist(aas))
  
})

for(i in 1:4){plot(graphs_min[[i]],main = abs[i])}


### select edges so that graphs are connected with minimal number of edges ----

disttemp = lapply(distances,\(df){
  df[df$Distance<=6,] 
  
})

distscan = lapply(distances,\(df){ # between 6 and 12
  df[df$Distance >6 & df$Distance <=12 ,]
  
})
distscan = lapply(distscan,\(df){ # between 6 and 12
  o = order(df$Distance) # put in ascending order
  df[o,]
  
})


gtemp = graphs_min # graphs to add edges to



# iterate through nodes over 6 A to include the necessary ones for a connected graph
for(i in 1:4){
  df = distscan[[i]]
  for(j in 1:nrow(df)){
    vs = V(gtemp[[i]])$name
    row = df[j,]
    if((row$Residue1 %in% vs == F) | (row$Residue2 %in% vs == F) | components(gtemp[[i]])$membership[row$Residue1] != components(gtemp[[i]])$membership[row$Residue2] ){
      d = row$Distance
      # if the chosen edge connects 2 disconnected components, include it with a gap node in the middle
      print(d)
      if(d%/%6 == 1){ # 1 gap node
        xnum = nrow(disttemp[[i]])
        n1 = paste0("XXX",xnum) # "gap" node 
        disttemp[[i]] = rbind(disttemp[[i]],c(row$Residue1,n1,d))
        disttemp[[i]] = rbind(disttemp[[i]],c(n1,row$Residue2,d))
        
      }else if (row$Distance%/%6 == 2){ # 2 gap nodes
        xnum = nrow(disttemp[[i]])
        n1 = paste0("XXX",xnum)
        n2 = paste0("XXX",xnum+1)
        disttemp[[i]] = rbind(disttemp[[i]],c(row$Residue1,n1,d))
        disttemp[[i]] = rbind(disttemp[[i]],c(n1,n2,d))
        disttemp[[i]] = rbind(disttemp[[i]],c(n2,row$Residue2,d))
      }
      gtemp[[i]] = graph_from_data_frame(disttemp[[i]],directed = F)
    }
  }
  
}
for(i in 1:4){plot(gtemp[[i]],main = abs[i])}

g_full = gtemp

# make AA 1 letter code as attribute
g_full = lapply(g_full,\(g){
  aas = sapply(V(g)$name,\(nm){
    c3 = paste0(unlist(strsplit(nm,""))[1:3],collapse = ""); 
    names(AA_code[which(AA_code==c3)])
    
  } )
  set_vertex_attr(g,name = "AA", value = aas)
  
})
save(g_full,file = "epitope_graphs.RData")

library(future.apply)
library(vctrs)
max_scores_all = future_sapply(1:4,\(ab_i){
  paths = longpaths(g_full[[ab_i]])
  sapply(all_pssm_,\(mot){
    
    scores = unlist(sapply(paths,\(paths2) {
      paths2 = list_drop_empty(paths2)
      sapply(paths2,\(path){
        score(mot,path)
      })
    }))
    max(scores)
  })
})

# make plots ----

for(i in 1:4){
  l = (len[which(colcodes==i)])
  o = order(l)
  sc = max_scores_all[which(colcodes==i),i][o]
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
    y = (sc+1.5)
    text(x[qth], y[qth], labels = round(qt[qth],2), cex = 1.2)
    
  }
  
}


# check which parts of epitope are similar to the clusters

for(i in 1:4){
  
  pssm = all_pssm_[which(colcodes==i)]
  
  sc = max_scores_all[which(colcodes==i),i]
  control = max_scores_all[,i]
  e = ecdf(control)
  qt = e(sc)
  qth = which(qt>0.8)
  
  paths = longpaths(g_full[[i]])
  printpaths = unlist(sapply(paths,\(paths2) {
    paths2 = list_drop_empty(paths2)
    sapply(paths2,\(path){
      paste(path,collapse = "")
    })
  }))
  for(q in qth){
    mot = pssm[[q]]
    scores = unlist(sapply(paths,\(paths2) {
      paths2 = list_drop_empty(paths2)
      sapply(paths2,\(path){
        score(mot,path)
      })
    }))
    jj = which.max(scores)
    print(paste(abs[i],q,round(scores[jj],2), printpaths[[jj]]))
  
  }
  
  
}
