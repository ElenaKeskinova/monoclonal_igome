AA_code = sapply(AMINO_ACID_CODE,toupper)
AA_code = c(AA_code,"XXX")
names(AA_code)[27] = "-"
abs = c("Herceptin","21c","17b","b12")

epi_contacts = lapply(abs,\(ab){
  read.delim(file = paste0("dock_motifier_mimotopes/",ab,"_contacts_epi.txt"))
})
epi_contacts = lapply(epi_contacts,\(ab){
  ab$AB_contacts = lapply(ab$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ",")); aas})
  ab
})
names(epi_contacts) = abs

# data frame of distances
pdistances = lapply(abs,\(ab){ # distances between amino acids closer than 5 A to antibody, from PDB, calculated with Pymol
  df = read.table(paste0(ab,"_pardist.txt"),header = T, sep = ",")
  df$aan1 = sapply(df$Residue1, substr, 1,3)
  df$aan2 = sapply(df$Residue2, substr, 1,3)
  df = df[which(df$aan1 %in% AA_code & df$aan2 %in% AA_code),] # only valid amino acids
  df
})
names(pdistances) = abs
allaas = lapply(pdistances,\(df){
  alaa = unique(unlist(df[,1:2]))})

dist_mx = lapply(pdistances,\(df){
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

twodim_coord = lapply(dist_mx,\(dmat){
  c = cmdscale(dmat,k = 2)
  c
})

# plot paratope map
i = 4
cnt = which(allaas[[i]] %in% unlist(epi_contacts[[i]]$AB_contacts))
plot(twodim_coord[[i]], main = abs[[i]], cex = 1.6)
points(twodim_coord[[i]][cnt,],pch = 16, col = "red",cex = 2)
text(twodim_coord[[i]],labels = allaas[[i]], cex = 0.8)
