# Calculates the number of seqs of length l in shell i on the Hamming lattice
# using alphabet of size aa.

hamLatShell=function(l,i, aa=20){
  aa=aa-1
  aa^i*choose(l,i)
}