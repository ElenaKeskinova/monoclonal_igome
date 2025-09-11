# Generates, (plots), and returns a subgraph of the Hamming lattice of
# subsequences of length k, corresponding to one aa sequence s

pepprint=function(s,k=5, pl=T){
  require(stringdist)
  require(igraph)
  
  x=unlist(strsplit(s, split=""))
  n=nchar(s)
  ii=combn(n,k)
  subStrs=unique(apply(ii,2,\(ij) paste(x[ij], collapse="")))
  dii=stringdistmatrix(subStrs, method = "ham")
  dii=as.matrix(dii); colnames(dii)=subStrs; rownames(dii)=subStrs
  dii[dii>1]=0
  g=graph_from_adjacency_matrix(dii, mode="undirected", add.colnames = NULL)
  if (pl) plot(g)
  return(g)
}