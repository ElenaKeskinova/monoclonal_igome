require(RBGL)
require(Biostrings)
require(parallel)
require(pbapply)
require(future.apply)
require(stringdist)
require(igraph)
require(reshape2)

aa=AA_STANDARD
L5=allpp("",aa,5)

length(unique(L5))
px=sample(L5,1)
pxd=stringdistmatrix(px,L5, method = "ham")
sh=table(pxd); sh=sapply(seq_along(sh),\(i) sum(sh[1:i]))
x=1:(length(sh)-1); sh=sh[-1];sh=cbind(x,sh)
# Correlation dimension
lm(formula = log(sh[2:4, 2]) ~ log(sh[2:4, 1]))$coefficients[2]

# formula for the N strings in shell i (if l is the number of letters, here 5) 
# N==19^i*combinations(l,i)=19^i*l!/((l-i)!*i!) (hamLatShell()), 

g=pepprint("ACDEFGH",5)

d=distances(g)
td=apply(d,1,table)
tdcd=sapply(td,\(x){
  xn=as.numeric(names(x))
  y=sapply(seq_along(x),\(i) sum(x[1:i]))
  lm(log(y[2:3])~log(xn[2:3]))$coefficients[2]
})
mean(tdcd)
ttd=table(d)
xn=as.numeric(names(ttd))
y=sapply(seq_along(ttd),\(i) sum(ttd[1:i]))
lm(log(y[2:3])~log(xn[2:3]))$coefficients[2]


z=allpp("",c(1,2,3,4,5,6,7),7)
zz=pbsapply(z,\(zi){
  zi=unlist(strsplit(zi,split=""))
  zu=unique(zi)
  zx=seq_along(zu)
  names(zx)=zu
  paste(zx[zi], collapse="")
})
names(zz)=NULL
zz=unique(zz)

gzz=pblapply(zz, \(zi){
  pepprint(zi,5)
})

isoj=pbsapply(gzz,\(gi){
  (sapply(gzz,\(gj){
    is_isomorphic_to(gj,gi)
  }))
})
giso=graph_from_adjacency_matrix(isoj*1, mode="undirected", diag = F)
jiso=components(giso)$membership
footprintClass=aggregate(zz, by=list(jiso), c)
for (ci in 1:114){
  j=grep(ci, components(giso)$membership)[1]
  plot(gzz[[j]], main=paste(ci,j,sep="_"))
}

xyneat=t(sapply(1:100,\(i) c(i,sum(1:i)*4)))
lm(log(xyneat[30:100,2])~log(xyneat[30:100,1]))$coefficients[2]
plot(log(xyneat[,1]),log(xyneat[,2]), ty="l")




iplan=sapply(gzz,\(gi) boyerMyrvoldPlanarityTest(as_graphnel(gi)))

g0=pepprint("ACDEFGH",5)
g1=pepprint("AIKEFGH",5)
g2=pepprint("ACYEFWH",5)
g3=pepprint("LCDEFGM",5)
g4=pepprint("ACDEFGM",5)
nms=lapply(list(g0,g1,g2,g3,g4), \(gi) vertex_attr(gi)$name); names(nms)=seq_along(nms)
nms=melt(nms);nms=aggregate(as.numeric(nms$L1), by=list(nms$value), \(x) if(length(x)>1) return(0) else return(x))
x=nms[,2];names(x)=nms[,1];nms=x
g0_4=graph.union(g0,g1, byname = T)
g0_4=graph.union(g0_4,g2, byname = T)
g0_4=graph.union(g0_4,g3, byname = T)
g0_4=graph.union(g0_4,g4, byname = T)
g0_4=set_vertex_attr(g0_4,"G",value=as.numeric(nms[vertex_attr(g0_4)$name]))
plot(g0_4, vertex.color=vertex_attr(g0_4)$G, vertex.label="", vertex.size=3)


