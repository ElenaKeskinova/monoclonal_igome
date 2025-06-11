library(Biostrings)
# use all_pssm_ from the script create_motifs.R
# take files with results from Pymol for contact amino acids on antigen surfaces, obtained from structures 1yyl,3lqa,2ny7,1n8z from pdb
# make graphs from epitope amino acids, introduce gap nodes where distance is over 6 A , take all paths and find highest overlap score to each cluster


AA_code = sapply(AMINO_ACID_CODE,toupper)
AA_code = c(AA_code,"XXX")
names(AA_code)[27] = "-"


distances = lapply(abs,\(ab){ # distances between amino acids closer than 5 A to antibody, from PDB, calculated with Pymol
  read.table(paste0("mixed-7graphs/epi_",ab,".txt"),header = T)
})


# clean from trash ----
aanames1 = lapply(distances,\(df){
  name1 =sapply(df$Residue1,\(r1){
    paste0(unlist(strsplit(r1,""))[1:3],collapse = "")
  })
  name1
  
})
aanames2 = lapply(distances,\(df){
  name2 =sapply(df$Residue2,\(r1){
    paste0(unlist(strsplit(r1,""))[1:3],collapse = "")
  }) 
  name2
  
})

for(i in 1:4){
  distances[[i]] = distances[[i]][which(aanames1[[i]]!="HOH" & (aanames2[[i]]!="HOH") & aanames1[[i]]!="NAG" & (aanames2[[i]]!="NAG")),]

}

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
  graph_from_edgelist(edges, directed = F)
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
library(foreach)


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


# max scores all graphs vs all clusters ----


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


# functions to compute paths and scores ----

longpaths = function(g){ # graph
  vs = V(g)[which(degree(g)==1)]
  peps = sapply(vs,\(v){
    paths = all_simple_paths(g,from =v,to = vs)
    peps = (sapply(paths,\(path){
      p = V(g)[names(path)]$AA
      np = length(p)
      nx = length(which(p=="-"))
      if(np-nx>4 & nx<np/3 ){ 
        p
      }
      # returns all paths as arrays
      
    }))
  })
  
}

score = function(mot,epi){ # provide a pssm and an epitope as a character array
  n = length(epi)
  if(n<=ncol(mot)){
    max(sapply(1:(ncol(mot)-n+1),\(i){
      sum(diag(mot[epi,i:(i+n-1)]))
    }))
  } 
  else{
    max(sapply(1:(n - ncol(mot)+1),\(i){
      sum(diag(mot[epi[i:(i+ncol(mot)-1)],]))
    }))
  }
}


