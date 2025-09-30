# amino acid code
library(future.apply)
library(Biostrings)
library(igraph)
AA_code = sapply(AMINO_ACID_CODE,toupper)
AA_code = c(AA_code,"XXX")
names(AA_code)[27] = "-"

# load set of epitope aminoacids

files = list.files(path = "pdb_dist", pattern = "^[A-Za-z0-9]{4}\\.txt$", full.names = TRUE)

df_dist = lapply(files,\(f){
  df = read.table(f,sep = ",",header = T)
  df$aan1 = sapply(df$Residue1, substr, 1,3)
  df$aan2 = sapply(df$Residue2, substr, 1,3)
  df[which(df$aan1 %in% AA_code & df$aan2 %in% AA_code),]

})
names(df_dist) = sapply(files,\(f) substr(f,10,13))

# create distance matrix

dist_mats = lapply(df_dist,\(df){
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
dshort = lapply(df_dist,\(df){
  df$Distance<=d
  
})

graphs_min_t = lapply(dist_mats,\(mat){
  require(igraph)
  
  graph_from_adjacency_matrix(mat>0 & mat<=d,mode = "undirected")
})
names(graphs_min_t) = names(df_dist)
for(i in names(graphs_min_t)){plot(graphs_min_t[[i]],main = i)}


### select edges so that graphs are connected with minimal number of edges ----

# disttemp = lapply(df_dist,\(df){
#   df[df$Distance<=6,] 
#   
# })
# 
# distscan = lapply(df_dist,\(df){ # between 6 and 12
#   df[df$Distance >6 & df$Distance <=12 ,]
#   
# })
# distscan = lapply(distscan,\(df){ 
#   o = order(df$Distance) # put in ascending order
#   df[o,]
#   
# })
# 
# 
# gtemp = graphs_min # graphs to add edges to
# 
# 
# 
# # iterate through nodes over 6 A to include the necessary ones for a connected graph
# for(i in names(gtemp)){
#   df = distscan[[i]]
#   for(j in 1:nrow(df)){
#     vs = V(gtemp[[i]])$name
#     row = df[j,]
#       d = row$Distance
#       if(d%/%6 <= 2){ # 1 gap node
#         xnum = nrow(disttemp[[i]])
#         n1 = paste0("XXX",xnum) # "gap" node 
#         disttemp[[i]] = rbind(disttemp[[i]],c(row$Residue1,n1,d))
#         disttemp[[i]] = rbind(disttemp[[i]],c(n1,row$Residue2,d))
#         
#       }else if (row$Distance%/%6 > 2){ # 2 gap nodes
#         xnum = nrow(disttemp[[i]])
#         n1 = paste0("XXX",xnum)
#         n2 = paste0("XXX",xnum+1)
#         disttemp[[i]] = rbind(disttemp[[i]],c(row$Residue1,n1,d))
#         disttemp[[i]] = rbind(disttemp[[i]],c(n1,n2,d))
#         disttemp[[i]] = rbind(disttemp[[i]],c(n2,row$Residue2,d))
#       }
#       gtemp[[i]] = graph_from_data_frame(disttemp[[i]],directed = F)
#     
#   }
#   
# }
# 
# for(i in names(gtemp)){plot(gtemp[[i]],main = paste(i,"full"))}
# 
# g_full = gtemp
# # make AA 1 letter code as attribute
# g_full = lapply(g_full,\(gf){
#   aas = sapply(V(gf)$name,\(nm){
#     c3 = substr(nm,1,3) 
#     names(AA_code[which(AA_code==c3)])
#     
#   } )
#   set_vertex_attr(gf,name = "AA", value = unlist(aas))
#   
# })

# 
# ## generate all paths of each antibody----
# source("D:/Elena/ban/monoclonal_igome/paths_scores_functions.R")
# library(vctrs)
# allpaths_control = lapply(g_full,\(gf) {
#   longpaths(gf)
# })

# load mimotopes ----
files_mimo = list.files(path = "pdb_dist", pattern = "^[A-Za-z0-9]{4}_mimo\\.txt$", full.names = TRUE)
mimotest = lapply(files_mimo, readLines)
names(mimotest) = sapply(files_mimo, substr,10,13)
mimotest_nC = lapply(mimotest,\(set){
  sapply(set,\(p){
    np = nchar(p)
    pp = unlist(strsplit(p,""))
    if((pp[1] == "C") & (pp[np] == "C")){
      substr(p,2,(np-1))
    }
    else{
      p
    }
  })
})

mimotest_nC7 = lapply(mimotest,\(set){
  unlist(sapply(set,\(p){
    np = nchar(p)
    pp = unlist(strsplit(p,""))
    if(np>7){
      sapply(1:(np-6),\(i){
        paste0(pp[i:(i+6)],collapse = "")
        
      })
    }
    else{
      p
    }
  }))
})


# # scores between mimotopes and paths----
# 
# scores_control_pos = lapply(names(mimotest_nC7),\(pdbname){
#   sapply(allpaths_control[[pdbname]],\(path){
#     sapply(mimotest_nC7[[pdbname]],\(mimo){
#       score_pep(mimo,path,blosum)
#     })
#     
#   })
#   
# })
# scores_pos_max = sapply(scores_control_pos,\(scmat){
#   apply(scmat,1,max)
# })
# names(scores_pos_max) = names(mimotest_nC)
# hist(unlist(scores_pos_max))
# 
# ## negative control----
# 
# # probabiliteas of ab-ag interaction
# ag_aaprob = read.csv(file = "AgAbIFprobs.csv")
# aa_prob = ag_aaprob$p
# names(aa_prob) = names(AMINO_ACID_CODE[1:20])
# # rnd peptides
# n = 200
# rnd_p = sapply(1:n,\(i){
#   
#   paste0(sample(AA_STANDARD,7,prob = aa_prob,replace = T),collapse = "")
# })
# 
# # scores
# scores_control_neg = future_lapply(names(mimotest_nC),\(pdbname){
#   sapply(allpaths_control[[pdbname]],\(path){
#     sapply(rnd_p,\(mimo){
#       score_pep(mimo,path,blosum)
#     })
#     
#   })
#   
# })
# scores_neg_max = sapply(scores_control_neg,\(scmat){
#   apply(scmat,1,max)
# })
# hist(unlist(scores_neg_max))
# 

# scores by following the graph and matching aminoacids like the MimoTree algorithm----

## prepare graph ----
graphs_min_t = lapply(graphs_min_t,\(g){ # 1 letter code as attribute
  aas = sapply(V(g)$name,\(nm){
    c3 = substr(nm,1,3) 
    names(AA_code[which(AA_code==c3)])
    
  } )
  set_vertex_attr(g,name = "AA", value = unlist(aas))
  
})

# search all partial overlaps like in  Mimotree, 
# list names are indices of starting amino acid in the peptide
# names of path elements are their scores, paths of last amino acid have scores as names of the whole list element!
# to do:
# join partial paths to bigger ones
# calculate overall scores as sum of individual
# compare scores of test sets to random -> roc curve to get treshold, compare to scores with blosum, chemical properties and binding correlation
source("D:/Elena/ban/monoclonal_igome/recursive_path_search.R")
### test with one mimotope ----
pep = "YEFQHY"
g = graphs_min_t[[1]]
paths = pepmatch(pep,g,blosum)

paths2 = merged_paths(paths,pep,dist_mats[[1]])

scores1 = sapply(paths2,\(p) sum(as.numeric(names(p))))
hist(scores1)

## all positive control scores ----

plan(multisession(workers = 4))
allscores = future_lapply(names(mimotest_nC7),\(pdbname){
  scores = sapply(mimotest_nC7[[pdbname]],\(mimo){
    paths = pepmatch(pep = mimo,graphs_min_t[[pdbname]],blosum)
    paths = merged_paths(paths,mimo,dist_mats[[pdbname]])
    
    max(sapply(paths,\(p) sum(as.numeric(names(p)))), na.rm = T) # take max score of all paths for each mimotope
    
    })
  print(pdbname)
  scores
})
hist(unlist(allscores))
## negative controls with random peptides ----

# probabiliteas of ab-ag interaction
ag_aaprob = read.csv(file = "AgAbIFprobs.csv")
aa_prob = ag_aaprob$p
names(aa_prob) = names(AMINO_ACID_CODE[1:20])

# rnd peptides
n = 100
rnd_p = sapply(1:n,\(i){
   paste0(sample(AA_STANDARD,7,prob = aa_prob,replace = T),collapse = "")
})
# with library background frequencies, see librarypeps.R
n = 100
rnd_p = sapply(1:n,\(i){
  paste0(sample(AA_STANDARD,7,prob = freql[AA_STANDARD],replace = T),collapse = "")
})

# same 100 with all epitopes
plan(multisession(workers = 6))
neg_scores = future_lapply(names(mimotest_nC7),\(pdbname){
  scores = sapply(rnd_p,\(mimo){
    paths = pepmatch(pep = mimo,graphs_min_t[[pdbname]],blosum)
    paths = merged_paths(paths,mimo,dist_mats[[pdbname]])
    
    max(sapply(paths,\(p) sum(as.numeric(names(p)))), na.rm = T) # take max score of all paths for each mimotope
    
  })
  print(pdbname)
  scores
})
names(neg_scores) = names(mimotest_nC7)
hist(unlist(neg_scores))
# different for each
# check 
hist(unlist(neg_scores),xlim = c(0,max(unlist(allscores))),col = rgb(1,0,0,0.5),main = "negative and positive controls")
hist(unlist(allscores),col = rgb(0,1,0,0.5),add = T)

library(pROC)
labels = c(rep(0,length(unlist(neg_scores))),rep(1,length(unlist(allscores))))
scores = c(unlist(neg_scores),unlist(allscores))
roc_obj <- roc(labels, scores)
# Plot ROC curve
plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve")
opt <- coords(roc_obj, "best", best.method = "closest.topleft")
opt_cutoff = opt$threshold
opt


# save data for roc curve 
save(mimotest_nC7,graphs_min_t,dist_mats,rnd_p,file = "testset_epimatch.R")
