
# all mimototes vs all paths for each antibody
allpepsets_w = lapply(abs, \(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  load(file = paste0(path,ab,"_peps_w.RData"))
  pepsets
})


all_mim_vs_epi = future_lapply(1:4,\(i){
    sapply(allpepsets_w[[i]],\(mimosets){
      sapply(mimosets,\(mimo){
        paths = pepmatch(pep = mimo,graphs_min[[i]],zcorr)
        paths = merged_paths(paths,mimo,dist_mx[[i]])
      
        max(sapply(paths,\(p) sum(as.numeric(names(p)))), na.rm = T) # take max score of all paths for each mimotope
      
      })
      
    })
})

