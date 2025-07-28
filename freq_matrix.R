# construct position probability matrix from alignment
freq_matrix = function(pep_arr, all_aa,ps_c = 0){
  require(Biostrings)
  m1 = consensusMatrix(pep_arr )
  mfull = matrix(0, nrow = length(all_aa), ncol = ncol(m1), dimnames = list(all_aa,colnames(m1)))
  
  if(rownames(m1)[1]=="-"){
    mfull[rownames(m1)[-1],] = m1[rownames(m1)[-1],] # for alignment, not peps
  
  }
  else{
    mfull[rownames(m1),] = m1[rownames(m1),] # for alignment, not peps
  }
  
  mfull = (mfull + ps_c/20)/(length(pep_arr)+ps_c)
  
  return(mfull)
}


# construct a position specific scoring matrix from alignment
pssm = function(peps,bkg_p = rep(1/20,20),aa = sort(Biostrings::AA_STANDARD),ps_c = 1){
  m = freq_matrix(peps,aa,ps_c)
  pssm = apply(m,2,\(cm)log2(cm/bkg_p))
  pssm
}
