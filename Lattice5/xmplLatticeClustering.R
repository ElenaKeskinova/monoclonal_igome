
require(Biostrings)
require(parallel)
require(pbapply)
require(future.apply)
require(stringdist)
require(igraph)
require(reshape2)
require(matrixStats)
require(Biostrings)
require(uwot)
require(rgl)
require(scales)
require(ggseqlogo)
pslg=pseudo_log_trans(sigm=1)
aa=AA_STANDARD

herpp=read.table("/home/anastas/Documents/EelenaK/monoclonal_igome/mixed-7graphs/Herceptin/Herceptin_allp810_noC7-fr.txt")
pp17b=read.table("/home/anastas/Documents/EelenaK/monoclonal_igome/mixed-7graphs/17b/17b_allp810_noC7.txt")
b12pp=read.table("/home/anastas/Documents/EelenaK/monoclonal_igome/mixed-7graphs/b12/b12_allp810_noC7.txt")
rndpp=sapply(1:8000,\(i) paste(sample(aa,7), collapse="", sep=""))
herpp=herpp[,1];pp17b=pp17b[,1];b12pp=b12pp[,1]
peps=c(herpp[1:2000],pp17b[1:2000], b12pp[1:2000], rndpp[1:2000])
pepsource=c(rep("her",2000),rep("17b",2000),rep("b12",2000),rep("rnd",2000))

ppMaps=pblapply(peps,\(p) {
  x=unlist(strsplit(p, split=""))
  n=nchar(p)
  ii=combn(n,5)
  unique(apply(ii,2,\(ij) paste(x[ij], collapse="")))
})

proct=proc.time()
n=detectCores()%/%3
plan(multisession, workers=n)
ppHausdorf=future_sapply(ppMaps,\(P){
    sapply(ppMaps, \(Q){
        M=stringdistmatrix(P, Q, method = "ham", nthread=3)
        mean(c(rowMins(M), colMins(M)))
      })
})
plan(sequential)
print(proc.time()-proct)

emb=cmdscale(scale(ppHausdorf),k=1400, eig = T, list.=T)
ume=umap(pslg$transform(emb$points), n_neighbors = 2000, n_components = 3, verbose=T, n_epochs=1000)
umed=umap(ppHausdorf, n_neighbors = 2000, n_components = 3, verbose=T, n_epochs=1000)

plot3d(ume, col=as.numeric(as.factor(pepsource)), size=3)
plot3d(emb$points[,1:3], col=as.numeric(as.factor(pepsource)), size=5)

rndppcoded=t(pbsapply(rndpp,\(p){
  x=unlist(strsplit(p,split=""))
  sapply(x,\(ci) which(aa==ci))  
}))

umeham=umap(rndppcoded, n_neighbors=50, n_components = 3, verbose=T, n_epochs=1000, metric="ham")
plot3d(umeham, size=3)

pepcoded=t(pbsapply(peps,\(p){
  x=unlist(strsplit(p,split=""))
  sapply(x,\(ci) which(aa==ci))  
}))

umehamall=umap(pepcoded, n_neighbors=50, n_components = 3, verbose=T, n_epochs=1000, metric="ham")
plot3d(umehamall, col=as.numeric(as.factor(pepsource)),size=3)

pepsL5coded=pblapply(peps,\(p){
  x=unlist(strsplit(p,split=""))
  n=nchar(p)
  ii=combn(n,5)
  l=unique(apply(ii,2,\(ij) paste(x[ij], collapse="")))
  t(sapply(l,\(p5){
    x=unlist(strsplit(p5,split=""))
    1*do.call(c,lapply(x,\(ai) aa==ai)) 
  }))
})

colL5coded=unlist(sapply(seq_along(pepsL5coded), \(i){
    l=pepsL5coded[[i]]
    rep(pepsource[i],nrow(l))
}))

pepsL5coded=do.call(rbind,pepsL5coded)
pepsL5coded=as.numeric(pepsL5coded)

umehall5=umap(pepsL5coded, n_neighbors = 50, n_components = 3, verbose=T, n_epochs=1000, metric="ham")

plot3d(umehall5, col=as.numeric(as.factor(colL5coded)), size=1)


# totally random set -----------------------------------------------------------


pprndMaps=pblapply(rndpp,\(p) {
  x=unlist(strsplit(p, split=""))
  n=nchar(p)
  ii=combn(n,5)
  unique(apply(ii,2,\(ij) paste(x[ij], collapse="")))
})

proct=proc.time()
n=detectCores()%/%3
plan(multisession, workers=n)
pprndHausdorf=future_sapply(pprndMaps,\(P){
  sapply(pprndMaps, \(Q){
    M=stringdistmatrix(P, Q, method = "ham", nthread=3)
    mean(c(rowMins(M), colMins(M)))
  })
})
plan(sequential)
print(proc.time()-proct)

emb=cmdscale(scale(ppHausdorf),k=1400, eig = T, list.=T)

umehamalL5=umap(pepcoded, n_neighbors=50, n_components = 3, verbose=T, n_epochs=1000, metric="ham")
plot3d(umehamall, col=as.numeric(as.factor(pepsource)),size=3)

# lattice-based graph ---------------------------------------------------------


