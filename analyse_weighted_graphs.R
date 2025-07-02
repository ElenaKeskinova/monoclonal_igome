abs = c("Herceptin","21c","17b","b12")
# create weighted graphs
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

# plot 3d graphics of clusters
ab = abs[1]
path = paste0("mixed-7graphs/",ab,"/")
load(paste0(path,ab,"big7-lw.RData"))
load(file = paste0(path,ab,"_",bc,"_","dbscan.RData"))
load(file = paste0(path,ab,"_",bc,"_","speccoord.RData"))

# make pepsets 


# make motifs


# compare to clusters without weights


# 