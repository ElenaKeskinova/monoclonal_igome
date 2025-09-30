# print a logo of a chosen set of peptides
printlogo = function(peps, nm = ""){
  require(ggseqlogo)
  require(ggplot2)
  require(Biostrings)
  require(msa)
  
  l=msaClustalW(AAStringSet(peps), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l@unmasked),1, paste,collapse="")
  
  m = consensusMatrix(l, as.prob=T)
  print(ggseqlogo(unlist(m))+
          ggtitle(nm))
  
}

logoal = function(peps, nm = ""){ # for aligned peptides
  require(ggseqlogo)
  require(ggplot2)
  require(Biostrings)
  
  m = consensusMatrix(peps, as.prob=T)
  print(ggseqlogo(unlist(m))+
          ggtitle(nm))
  
}