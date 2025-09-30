library(future.apply)
# functions---- 

get_contacts = function(path, pattern){
  dockfiles = list.files(path = path,pattern = "^model_cluster[0-9]+_[0-9]+\\.pdb_contacts\\.txt$")
  clname = sapply(dockfiles, \(name) unlist(strsplit(name,"_"))[2])
  uninames = sapply(dockfiles, \(name) paste0(unlist(strsplit(name,"_"))[2:3],collapse = "_"))
  allcontacts = lapply(dockfiles,\(f){
    df = read.delim(file = paste0(path,"/",f))
    df$AB_contacts = lapply(df$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ",")); sapply(aas, \(aa) substr(aa,1,nchar(aa)-1))})
    df$aaID = sapply(df$AG_contact,\(aa){as.numeric(substr(aa,4,nchar(aa)-1))})
    df
  })
  names(allcontacts) = uninames
  allcontacts
}

contact_matrix = function(clnames, epi_c, mimocontacts, mimo){ # absolute number of common contacts
  cluni = unique(clnames)
  nr = nrow(epi_c)
  
  matches = lapply(cluni,\(cl){
    nums = which(clnames ==cl)
    mimo_c = mimocontacts[nums]
    con_mat = matrix(0, nrow = nr, ncol = nchar(mimo) )
    rownames(con_mat) = epi_c$AG_contact
    colnames(con_mat) = unlist(strsplit(mimo,""))
    for(df in mimo_c){
      for(i in 1:nr){
        for(j in 1:nchar(mimo)){
          con_mat[i,j] = con_mat[i,j] + length(intersect(unlist(epi_c[i,"AB_contacts"]),unlist(df[which(df$aaID == j),"AB_contacts"])))
        }
        
      }
    }
    con_mat
    
  })
  names(matches) = cluni
  matches
}

contact_matrix_ = function(clnames, epi_c, mimocontacts, mimo){ # percentage of overlapping 
  cluni = unique(clnames)
  nr = nrow(epi_c)
  
  matches = lapply(cluni,\(cl){ # for each cluster
    nums = which(clnames ==cl)
    mimo_c = mimocontacts[nums]
    con_mat = matrix(0, nrow = nr, ncol = nchar(mimo) )
    
    
    con_mats = lapply( mimo_c,\(df){ # iterate over models in a cluster
      
      mat = sapply(1:nr,\(i){
        nr_epi = (unlist(epi_c[i,"AB_contacts"])) # contacts in epitope
        
        sapply(1:nchar(mimo),\(j){
          nr_all = length(unique(nr_epi,unlist(df[which(df$aaID == j),"AB_contacts"])))
          nr_common = length(intersect(unlist(epi_c[i,"AB_contacts"]),unlist(df[which(df$aaID == j),"AB_contacts"])))
          nr_common/nr_all
        })
        
      })
      mat = t(mat)
      rownames(mat) = epi_c$AG_contact
      colnames(mat) = unlist(strsplit(mimo,""))
      mat
    })
    Reduce("+", con_mats)/length(con_mats) # mean over all matrices for 1 cluster
    
    
  })
  names(matches) = cluni
  matches
}

###############################################
#--------------------------------------------------------------------
t_abs = c("1bj1","1iqd","1yy9","2adf")
testmimo = read.delim(file = "pep_to_dock.txt")


# all contacts between epitope and paratope
epi_contacts = lapply(t_abs,\(ab){
  read.delim(file = paste0("pepdock/",ab,"_contacts_epi.txt"))
})
epi_contacts = lapply(epi_contacts,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ",")); sapply(aas, \(aa) substr(aa,1,nchar(aa)-1))})
  ab
})
names(epi_contacts) = t_abs

gmin_test = graphs_min_t[t_abs]
gmin_test = lapply(t_abs,\(ab) subgraph(gmin_test[[ab]],which(V(gmin_test[[ab]])$name %in% epi_contacts[[ab]]$AG_contact)))

# extract contacts for all antibodies and mimotoopes
paths = c(sapply(t_abs,\(nm){paste0("pepdock/",nm,"_",c(1,2),"_docked")}))
contacts_df = sapply(paths, get_contacts)
names(contacts_df) = t_abs
mimos = as.vector(t(testmimo[,2:3]))
clnames_all = lapply(paths,\(path){
  dockfiles = list.files(path = path,pattern = "^model_cluster[0-9]+_[0-9]+\\.pdb_contacts\\.txt$")
  sapply(dockfiles, \(name) unlist(strsplit(name,"_"))[2]) }) 

plan(multisession(workers = 4))
overlap_matrices = future_lapply(seq_along(mimos),\(i){ # iterate over mimotopes
  clnames = clnames_all[[i]]
  epi_c = epi_contacts[[ceiling(i/2)]] # 2 mimotopes per antibody
  mimo = mimos[i]
  mimocont = contacts_df[[i]]
  contact_matrix_(clnames, epi_c, mimocont, mimo)
}) 
# put names
names(overlap_matrices) = sapply(1:8,\(i) paste0(t_abs[ceiling(i/2)],"_",mimos[i]))


# total number of common contacts
plan(multisession(workers = 4))
overlap_matrices_t = future_lapply(seq_along(mimos),\(i){ # iterate over mimotopes
  clnames = clnames_all[[i]]
  epi_c = epi_contacts[[ceiling(i/2)]] # 2 mimotopes per antibody
  mimo = mimos[i]
  mimocont = contacts_df[[i]]
  contact_matrix(clnames, epi_c, mimocont, mimo)
}) 
# put names
names(overlap_matrices_t) = sapply(1:8,\(i) paste0(t_abs[ceiling(i/2)],"_",mimos[i]))

# make matrices with strength and mean position of contact

scales = lapply(overlap_matrices,\(model){
  lapply(model,\(clusterm){
    intensity = rowSums(clusterm)
    position = sapply(1:nrow(clusterm),\(j) sum(clusterm[j,]*(1:ncol(clusterm)))/sum(clusterm[j,]))
    mat = cbind(intensity,position)
    rownames(mat) = rownames(clusterm)
    mat
  })
})

scales_t = lapply(overlap_matrices_t,\(model){
  lapply(model,\(clusterm){
    intensity = rowSums(clusterm)
    position = sapply(1:nrow(clusterm),\(j) sum(clusterm[j,]*(1:ncol(clusterm)))/sum(clusterm[j,]))
    
    mat = cbind(intensity,position)
    rownames(mat) = rownames(clusterm)
    mat
  })
})



# # add vertex attribute number of matches
# V(gmin_test[[1]])$nb_matches = 0
# for(aai in rownames(matches[[1]])){
#   V(gmin_test[[1]])[aai]$nb_matches=V(gmin_test[[1]])[aai]$nb_matches+rowSums(matches[[3]])[aai]
# }
# plot(gmin_test[[1]],vertex.size =  V(gmin_test[[1]])$nb_matches+1)

# project graphs in one dimension
onedim_scale = lapply(dist_mats[t_abs],\(dmat){
  
  s = cmdscale(dmat,k=1)
  names(s) = rownames(dmat)
  s = sort(s)
  names(s)
})

onedim_scale = lapply(t_abs,\(i){
  onedim_scale[[i]][which(onedim_scale[[i]] %in% epi_contacts[[i]]$AG_contact)]
})
names(onedim_scale) = t_abs

save(onedim_scale,matches,file = "dockresults_match_g1mimo1.RData")

# project graphs in two dimensions

twodim_scale = lapply(dist_mats[t_abs],\(dmat){
  
  s = cmdscale(dmat,k=2)
  rownames(s) = rownames(dmat)
  s
})
twodim_scale_cl = lapply(t_abs,\(i){
  twodim_scale[[i]][which(rownames(twodim_scale[[i]]) %in% epi_contacts[[i]]$AG_contact),]
})

cpl = colorRampPalette(c("blue","red")) # colors for contacts from front to back
for(i in 1:4){
  for(j in 1:2){
    nmb = (i-1)*2+j
    mimo = mimos[nmb]
    clnames = unique(clnames_all[[nmb]])
    pl = twodim_scale_cl[[i]]
    
    for(cl in clnames){
      sc = scales_t[[nmb]][[cl]]
      sc = sc[rownames(pl),]
      plot(pl,pch = 16, cex = log2(sc[,1]+1), col = cpl(nchar(mimo))[sc[,2]], xlim = range(pl[,1]) + c(-3,3),ylim = range(pl[,2]) + c(-3,3),main = paste(t_abs[i],mimo,cl))
      #plot(pl,pch = 16, cex = (sc[,1])*10, col = cpl(nchar(mimo))[sc[,2]], xlim = range(pl[,1]) + c(-3,3),ylim = range(pl[,2]) + c(-3,3),main = paste(t_abs[i],mimo,cl))
      
      text(pl+1, labels = rownames(pl))
    }
    
  }
  
}

# compare good scores from matrices with amino acid similarity matrices
library(reshape)
library(Biostrings)
AAmats = lapply(overlap_matrices, \(cls) lapply(cls,\(mat) melt(mat)) )
AAmats = (unlist(AAmats,recursive = F))
AAmats = do.call(rbind,AAmats)
AAmats$AA1 = sapply(AAmats$X1,\(aanm) {
  trilet = substr(aanm,1,3)
  names(AA_code)[which(AA_code == trilet)]
})
AAmats$AA2 = as.character(AAmats$X2)
AAmats$Zscalecor = sapply(1:nrow(AAmats), \(i) zcorr[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$blosum = sapply(1:nrow(AAmats), \(i) blosum[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$interface = sapply(1:nrow(AAmats), \(i) corrmat[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$Z1 = sapply(1:nrow(AAmats), \(i) Z1m[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$Z2 = sapply(1:nrow(AAmats), \(i) Z2m[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$Z3 = sapply(1:nrow(AAmats), \(i) Z3m[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$Z4 = sapply(1:nrow(AAmats), \(i) Z4m[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmats$Z5 = sapply(1:nrow(AAmats), \(i) Z5m[AAmats[i,"AA1"],AAmats[i,"AA2"]])
AAmat = AAmats[,c("AA1","AA2","value","Zscalecor","blosum","interface","Z1","Z2","Z3","Z4","Z5")] # reorder

AAmats_good = AAmat[which(AAmat$value>0.05),]
plot(AAmats_good$value,AAmats_good$Zscalecor)
plot(AAmats_good$value,AAmats_good$blosum)
plot(AAmats_good$value,AAmats_good$interface)
plot(AAmats_good$value,AAmats_good$Z3)


AAmats_good = AAmat[which(AAmat$value>10),]
model = lm(value~Zscalecor+blosum+interface+Z1+Z2+Z3+Z4+Z5,AAmats_good)
summary(model)
save(AAmat, file = "AAmat_corr.RData")
load(file = "AAmat_corr.RData")

# Initialize square matrix with 0
mat <- matrix(0, nrow = 20, ncol = 20,
              dimnames = list(AA_STANDARD,AA_STANDARD))

# Fill matrix with mean similarities
for (a1 in AA_STANDARD) {
  for(a2 in AA_STANDARD){
    rows <- which(AAmat$AA1 == a1 & AAmat$AA2 == a2 ) 
    
    mat[a1, a2] <- mean(AAmat[rows,"value"], na.rm = TRUE)

  }
  
}

mat
mat[which(is.na(mat))] = 0
heatmap(mat,scale = "none")

save(gmin_test, AAmat, twodim_scale_cl, overlap_matrices, scales, dist_mats, t_abs, file = "pepdock/docking_contacts_and_graphs.RData")
