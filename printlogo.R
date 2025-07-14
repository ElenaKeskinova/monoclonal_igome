# print a logo of a chosen set of peptides
printlogo = function(peps){
  require(ggseqlogo)
  require(ggplot2)
  require(Biostrings)
  
  l=msaClustalW(AAStringSet(peps), gapOpening = 2, gapExtension = 1, maxiters=1000, substitutionMatrix = "blosum")
  l= apply(as.matrix(l),1, paste,collapse="")
  
  m = consensusMatrix(l, as.prob=T)
  print(ggseqlogo(unlist(m)))
  
}
