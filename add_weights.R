# put weights on the edges of line graphs, according to the copy number of the respective sequences
add_weights = function(ab){
  require(future.apply)
  path = paste0("mixed-7graphs/",ab,"/")
  
  load(paste0(path,ab,"big7or.RData"))
  load(paste0(path,ab,"big7-l.RData"))
  load(file = paste0(path,ab,"_all_lcs.RData"))
  
  v_freqs = V(G)$Freq
  names(v_freqs) = V(G)$name
  
  lcs = V(lg)$name
  plan(multisession)
  weights = future_sapply(E(lg),\(e){
    ee = ends(lg,e)
    ee1 = ee[1]
    ee2 = ee[2]
    
    seqs1 = unlist(edg_lcs[which(edg_lcs[,3]==ee1),1:2])
    seqs2 = unlist(edg_lcs[which(edg_lcs[,3]==ee2),1:2])
    common = intersect(seqs1,seqs2)
    sum(v_freqs[common])
  })
  set_edge_attr(lg,name = "weight",value = weights)
}

weights = function(lg,edg_lcs,v_freqs){
  future_sapply(E(lg),\(e){
    ee = ends(lg,e)
    ee1 = ee[1]
    ee2 = ee[2]
    
    seqs1 = unlist(edg_lcs[which(edg_lcs[,3]==ee1),1:2])
    seqs2 = unlist(edg_lcs[which(edg_lcs[,3]==ee2),1:2])
    common = intersect(seqs1,seqs2)
    sum(v_freqs[common])
  })
}