# functions---- 

get_contacts = function(path, pattern= "^model_cluster[0-9]+_[0-9]+\\.pdb_contacts\\.txt$"){
  dockfiles = list.files(path = path,pattern )
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

read_cont = function(file){ # read file with contacts to a data frame
  df = read.delim(file)
  df$AB_contacts = lapply(df$AB_contacts,\(c_list){aas = unlist(strsplit(c_list, split = ",")); sapply(aas, \(aa) substr(aa,1,nchar(aa)-1))})
  df$aaID = sapply(df$AG_contact,\(aa){as.numeric(substr(aa,4,nchar(aa)-1))})
  df
}

cont_mat = function(epi_c, mimo_c, mimo){
  # 1 data frame with epi contacts, a list with data frames of mimotope contacts, and a string of the mimotope
  nr = nrow(epi_c)
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
