# Generates all aa seqs of a given length n. A is an empty set initially. 


allpp=function(A,aa,n){
  require(future.apply)
  plan(multisession, workers=n)
  A=c(future_sapply(A, \(Ai) sapply(aa,\(ai) paste(c(Ai,ai), collapse=""))))
  plan(sequential)
  if (nchar(A[1])<n) allpp(A,aa,n) else return(A)
}