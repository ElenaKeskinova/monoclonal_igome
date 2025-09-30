# compare kld to kld from the random clusters, calculate z-scores, make graphics


# get random results
load(file="rnd_len_kld_2.RData")
load( "rnd_pepsets_bkgfreqs_2.RData")

# calculate z scores from distribution of the kld of random clusters----
# z_sc = lapply(1:4,\(i){
#   data = data.frame(x=log10(rnd_results[[i]][,3]),y=log10(rnd_results[[i]][,1]))
#   real_l = loglengths[which(colcodes==i)]
#   names(real_l) = names(lengths[which(colcodes==i)])
#   real_k = log10(kl[which(colcodes==i)])
#   names(real_k)= names(kl[which(colcodes==i)])
#   model=kreg(data$x,data$y,grid = data$x)
#   pred =kreg(data$x,data$y,grid = real_l,bandwidth = model$bandwidth)
#   ii = order(data$x)
#   sq_res = (model$y-data$y[ii])**2
#   #hist(sq_res)
#   mod_res = kreg(model$x,sq_res,grid = data$x)#,bandwidth = 11.3
#   
#   pred_res = kreg(model$x,sq_res,grid = loglengths[which(colcodes==i)],bandwidth = mod_res$bandwidth) # var predictions for real data
#   sd = sqrt(pred_res$y)
#   
#   ij = order(real_l)
#   
#   z = as.vector((real_k[ij]-pred$y)/sd)
#   big = which(is.na(pred$y)) # out of the range of the randomly generated clusters
#   l = length(model$x)
#   
#  
#   m = which.max(pred$x[-big])
#   z[big] = as.vector((real_k[ij][big]-pred$y[-big][m])/sd[-big][m])
#   names(z) = names(sort(real_l))
#   z
#   
# })
library(scales)
library(gplm)
data = data.frame(log10(rnd_len_kld))
real_l = log10(lengths)
names(real_l) = names(lengths)
real_k = log10(all_kld)
names(real_k)= names(all_kld)
model=kreg(data$l,data$kl,grid = data$l)#Scott's rule of thumb for bandwidth 
pred =kreg(data$l,data$kl,grid = real_l,bandwidth = model$bandwidth)
ii = order(data$l)
sq_res = (model$y-data$kl[ii])**2

# calculate variance as prediction of residuals:
mod_res = kreg(model$x,sq_res,grid = data$l)
pred_res = kreg(model$x,sq_res,grid = real_l,bandwidth = mod_res$bandwidth) # variance predictions for real data
sd = sqrt(pred_res$y)

ij = order(real_l)
z_sc = as.vector((real_k[ij]-pred$y)/sd) #ordered
names(z_sc) = names(real_l[ij])
z_sc = z_sc[names(real_l)]

# p- values----
allp = sapply((z_sc),\(z){
  pnorm(z,lower.tail = F)
})
allp = p.adjust(allp,method = "BH") # ordered



# make graphics----
allnames = c(names(unlist(logpepsets,recursive = F)),"bkg")
### mapping distances

load(file = "4st_generation_library/lib_freqs.RData")
bgmot = freqM7


allp_mb = append(allm_wp2,bkg_m)
m_comp_b = compare_motifs(allp_mb,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
plot2d = cmdscale(m_comp_b,2)
rownames(plot2d) = c(allnames,"bkg")

sign_p  = allp<0.01
nn = nrow(plot2d)
# each antibody
for(i in 1:4){
  plot(plot2d[nn,1],plot2d[nn,2],cex = log((1000))/2,pch = 18,xlim =range( plot2d[,1]),ylim = range(plot2d[,2]),main = paste("clusters of",abs[i]),xlab = "dimension 1",ylab = "dimension 2" )
  
  num = which(numcodes == i)
  o = order(lengths[num])
  
  p = sign_p[num]
  
  #significant
  points(plot2d[num,][o,][p,],cex = log(sort(lengths[num][p]))/2,pch = 16,col = adjustcolor("orange", alpha.f = 0.7))
  #non-significant
  points(plot2d[num,][o,],cex = log(sort(lengths[num]))/2,col = adjustcolor("orange", alpha.f = 0.7))
}

# part of the poster:
# all anti gp120 antibodies together
colors = c("red","green","blue","orange")
plot(plot2d[nn,1],plot2d[nn,2],cex = log((1000))/2,pch = 18,xlim =range( plot2d[,1]),ylim = range(plot2d[,2]),main = paste("Clusters of anti-gp120\nantibodies"),xlab = "dimension 1",ylab = "dimension 2" )

for(i in 2:4){
  
  num = which(numcodes == i)
  o = order(lengths[num])
  
  p = sign_p[num]
  
  #significant
  points(plot2d[num,][o,][p,],cex = log(sort(lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot2d[num,][o,],cex = log(sort(lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
}
legend("topleft",legend = abs[2:4],col = colors[2:4],pt.cex = 1.6,cex = 1.3,pch = 16,title = "Antibodies",border = F)
legend("bottomleft",x.intersp = 2,pt.cex = log((1000))/4,cex = 1.3,pch = 18,legend = "background")



# only Herceptin
plot(plot2d[nn,1],plot2d[nn,2],cex = log((1000))/2,pch = 18,xlim =range( plot2d[,1]),ylim = range(plot2d[,2]),main = paste("Clusters of Herceptin"),xlab = "dimension 1",ylab = "dimension 2" )

for(i in 1:1){
  
  num = which(numcodes == i)
  o = order(lengths[num])
  
  p = sign_p[num]
  
  #significant
  points(plot2d[num,][o,][p,],cex = log(sort(lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot2d[num,][o,],cex = log(sort(lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
}
#legend("topleft",legend = abs[1],col = colors[1],pt.cex = 1.5,pch = 16,title = "Antibodies",border = F)
legend("bottomleft",x.intersp = 2,pt.cex = log((1000))/4,cex = 1.3,pch = 18,legend = "background")
points(plot2d["Herceptin_27",1],plot2d["Herceptin_27",2],pch = 4)


#only legend
plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)

legend("left",legend = c("p>0.01","p<0.01"),col =adjustcolor(colors[1], alpha.f = 0.7),pt.cex = 1.3,pch = c(1,16),title = "p-value",border = F)
legend("center",legend = c("10","50","100","500"),col = adjustcolor(colors[1], alpha.f = 0.7),pch = 16,pt.cex =log(c(10,50,100,500))/2,cex = 1,title = "Number of peptides \nper cluster",x.intersp = 1,y.intersp = 1.2,text.width = 0.15 )

plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
legend("center",legend = abs,col = colors,pt.cex = 1.6,cex = 1.3,pch = 16,title = "Antibodies",border = F)

# calculate Bhattacharyya coefficient and distance ----
allnames = c(names(unlist(logpepsets,recursive = F)),"bkg")
all_motifs = sapply(unlist(all_ppm,recursive = F), create_motif)
allp_mb = append(all_motifs,bkg_m)
m_comp_b = compare_motifs(allp_mb,tryRC = FALSE,min.mean.ic = 0.01,method = "BHAT")
m_comp_b = -log(m_comp_b)
colnames(m_comp_b) = allnames


# 3d pictures all abs together ----

allp_mb = append(allm_w,bkg_m)
m_comp_b = compare_motifs(allp_mb,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
plot_3d = cmdscale(m_comp_b,20, eig = T)
## eigv
sum = 0
eig = plot_3d$eig
eig = eig[eig>0]
eigvar = eig/sum(eig)
eigsum = sapply(1:length(eigvar),\(j) {
  sum <<- eigvar[j] + sum
  sum
  })
plot(eigsum)
#####

plot_3d = plot_3d$points
rownames(plot_3d) = c(names(unlist(logpepsets,recursive = F)),"bkg")

psl = transform_pseudo_log(sigma=0.01)
plot_3d = psl$transform(plot_3d)

sign_p  = allp<0.01
nn = nrow(plot_3d)

# all  antibodies together
colors = c("red","green","blue","orange")

d1 = 1
d2 = 3
plot(plot_3d[nn,d1],plot_3d[nn,d2],cex = log((1000))/2,pch = 18,xlim =range( plot_3d[,1]),ylim = range(plot_3d[,2]),main = paste("Clusters of all antibodies"),xlab = paste("dimension",d1),ylab = paste("dimension",d2) )

for(i in 1:4){
  
  num = which(numcodes == i)
  count = length(num)
  
  p = sign_p[allnames[num]]
  
  #significant
  points(plot_3d[num,c(d1,d2)][p,],cex = log((lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot_3d[num,c(d1,d2)],cex = log((lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
  
  #top 3 matches to epitope
  top = names((epi_p[[i]][count]))
  #top3 = names((epi_p[[i]][(count-2):(count)]))
  points(plot_3d[top,c(d1)],plot_3d[top,c(d2)],cex = log((lengths[top3]))/3,pch = 0,col = adjustcolor("black", alpha.f = 1))
}
legend("topright",legend = abs[1:4],col = colors[1:4],pt.cex = 1.6,cex = 1.3,pch = 16,title = "Antibodies",border = F)
legend("bottomright",x.intersp = 2,pt.cex = log((1000))/4,cex = 1.3,pch = 18,legend = "background")

plot3d(plot_3d,col = c(colors[numcodes],"black"),size = 5)

# project in a lower dimension ----
library(uwot)
p3d = umap(plot_3d)
p = sign_p
plot(p3d, col = c(colors[numcodes],"black"),cex = log(c(lengths,1000))/3)
points(p3d, col = adjustcolor(colors[numcodes[p]],0.3),cex = log(c(lengths[p]))/3,pch = 16)

d1 = 1
d2 = 2
plot(p3d[nn,d1],p3d[nn,d2],cex = log((1000))/2,pch = 18,xlim =range( p3d[,d1]) + c(0.1, 1),ylim = range(p3d[,d2])+ c(-0.5, 0.5),main = paste("Clusters of all antibodies"),xlab = paste("dimension",d1),ylab = paste("dimension",d2) )

for(i in 1:4){
  
  num = which(numcodes == i)
  count = length(num)
  
  p = sign_p[allnames[num]]
  
  #significant
  points(p3d[num,c(d1,d2)][p,],cex = log((lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(p3d[num,c(d1,d2)],cex = log((lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
  
  #top 3 matches to epitope
  top = names((epi_p[[i]][count]))
  #top3 = names((epi_p[[i]][(count-2):(count)]))
  points(p3d[top,c(d1)],p3d[top,c(d2)],cex = log((lengths[top3]))/3 + 1,pch = 0,col = adjustcolor("black", alpha.f = 1))
}
legend("right",legend = abs[1:4],col = colors[1:4],pt.cex = 1,cex = 1,pch = 16,title = "Antibodies",border = "white",bty = "n")
legend("bottomright",x.intersp = 2,pt.cex = log((1000))/4,cex = 1,pch = 18,legend = "background",border = "white",bty = "n")
legend("topright",legend = c("10","100","500"),col = adjustcolor(colors[1], alpha.f = 0.7),pch = 16,pt.cex =log(c(10,100,500))/3,cex = 1,title = "Number of \npeptides \nper cluster",x.intersp = 1,y.intersp = 1.2,text.width = 0.15,border = "white",bty = "n" )

#only legend
plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)

legend("left",legend = c("p>0.01","p<0.01"),col =adjustcolor(colors[1], alpha.f = 0.7),pt.cex = 1.3,pch = c(1,16),title = "p-value",bty = "n")
legend("center",legend = c("10","100","500"),col = adjustcolor(colors[1], alpha.f = 0.7),pch = 16,pt.cex =log(c(10,100,500))/3,cex = 1,title = "Number of peptides \nper cluster",x.intersp = 1,y.intersp = 1.2,text.width = 0.15 )
legend("right", pt.cex = log((1000))/4,cex = 1,pch = 0,legend = "closest to epitope" ,bty = "n")
#############3
rndpep = sapply(1:1000,\(i){
  paste(sample(sort(AA_STANDARD), 7,prob = rowSums(bgmot)/7, replace = TRUE), collapse = "")
})
rndpepm = sapply(rndpep,\(p) create_motif(p,alphabet = "AA"))

rndpepm = (future_sapply(rnd_pepsets,\(set){
  create_motif(set,alphabet = 'AA')
  }))
allp_rnd = append(allm_w,list(bkg_m,rndpepm))
allp_rnd = unlist(allp_rnd,recursive = F)
nnr =length(allp_rnd)
m_comp_all = compare_motifs(allp_rnd,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
cmd = cmdscale(m_comp_all,nnr-1)
rownames(cmd)[1:nn] = c(names(unlist(logpepsets,recursive = F)),"bkg")


library(FactoMineR)
pcar = PCA(cmd[(nn+1):nnr,],ncp=2)
colnames(cmd) = colnames(pcar$call$X)
pcall = predict(pcar,newdata = cmd)
pcpl = pcall$coord

# pcall = PCA(cmd,ncp=2)
# pcpl = pcall$ind$coord


plot(pcar$ind$coord[,1:2],cex = log((10))/2,pch = 0,main = paste("Clusters of all antibodies"),xlab = paste("dimension",d1),ylab = paste("dimension",d2),col = adjustcolor("grey", alpha.f = 0.3),xlim = range(pcar$ind$coord[,d1]),ylim = range(pcar$ind$coord[,d2]) )
points(pcpl[nn,d1],pcpl[nn,d2],cex = log(1000)/2,pch = 18,col = "black" )
for(i in 1:4){
  num = which(numcodes == i)
  p = sign_p[allnames[num]]
  #significant
  points(pcpl[num,c(d1,d2)][p,],cex = log((lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot_3d[num,c(d1,d2)],cex = log((lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
  
  }
legend("topleft",legend = abs[1:4],col = colors[1:4],pt.cex = 1.6,cex = 1.3,pch = 16,title = "Antibodies",border = F)
legend("bottomleft",x.intersp = 2,pt.cex = log((1000))/4,cex = 1.3,pch = 18,legend = "background")



allp_rnd_sep = lapply(1:4,\(i) unlist(list(allp_rnd[which(numcodes == i)],allp_rnd[94:nnr]),recursive = F))

cmd_all = lapply(allp_rnd_sep,\(mots)  cmdscale(compare_motifs(mots,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC"),2))

for(i in 1:4){
  b = table(numcodes)[i]+1
  plot_3d = cmd_all[[i]]
  num = which(numcodes == i)
  count = length(num)
  rownames(plot_3d)[1:count] = allnames[num]
  
  
  
  p = sign_p[allnames[num]]
  
  # random
  plot(plot_3d[(b+1):nrow(plot_3d),c(d1,d2)],cex = log((10))/2,pch = 0,xlab = paste("dimension",d1),ylab = paste("dimension",d2),col = adjustcolor("grey", alpha.f = 0.4),xlim =range( plot_3d[,1]),ylim = range(plot_3d[,2]),main = abs[i])
  # background
  points(plot_3d[b,d1],plot_3d[b,d2],cex = log((100))/2,pch = 18,col = "black" )

  #significant
  points(plot_3d[1:(b-1),c(d1,d2)][p,],cex = log((lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot_3d[1:(b-1),c(d1,d2)],cex = log((lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
  
  #top3 = names((epi_p[[i]][(count-2):(count)]))
  #points(plot_3d[top3,c(d1,d2)],cex = log((lengths[top3]))/3,pch = 8,col = adjustcolor("black", alpha.f = 1))
}

dist_all = sapply(1:4, \(i){
  coord = cmd_all[[i]][which(numcodes == i),]
  d = sqrt(rowSums(coord**2))
  names(d)=allnames[which(numcodes == i)]
  d
  
})

dist_sign = sapply(1:4, \(i){
  num = which(numcodes == i)
  count = length(num)
  d = dist_all[[i]]
  d[names(which(epi_p[[i]][(count-4):(count)]<0.1))]
  
})
dist_sign = list_drop_empty(dist_sign)

dist_unsign = sapply(1:4, \(i){
  num = which(numcodes == i)
  count = length(num)
  d = dist_all[[i]]
  d[names((epi_p[[i]][1:(count-5)]))]
  
})
boxplot(unlist(dist_unsign),unlist(dist_sign),varwidth = T,notch = T)

plot_3d = cmd
# random
  plot(plot_3d[95:nrow(plot_3d),c(d1,d2)],cex = log((10))/2,pch = 0,xlab = paste("dimension",d1),ylab = paste("dimension",d2),col = adjustcolor("grey", alpha.f = 0.4),xlim =range( plot_3d[,1]),ylim = range(plot_3d[,2]),main = abs[i])
  # background
  points(plot_3d[94,d1],plot_3d[94,d2],cex = log((100))/2,pch = 18,col = "black" )

for(i in 1:4){
  b = table(numcodes)[i]+1
  plot_3d = cmd_all[[i]]
  
  num = which(numcodes == i)
  count = length(num)
  
  p = sign_p[allnames[num]]
  
  
  #significant
  points(plot_3d[num,c(d1,d2)][p,],cex = log((lengths[num][p]))/3,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot_3d[num,c(d1,d2)],cex = log((lengths[num]))/3,col = adjustcolor(colors[i], alpha.f = 0.7))
  
}

  
  # graph from clusters and modularity ----
dist_matrix = compare_motifs(allm_w2,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
rownames(dist_matrix) = allnames[1:88]
colnames(dist_matrix) = allnames[1:88]

# select cutoff for graph:
hist(dist_matrix)
adj_matrix = dist_matrix<=1
g_mots = graph_from_adjacency_matrix(adj_matrix,mode = "undirected",diag = F)
edgedist = sapply(E(g_mots),\(e) {
  ee = ends(g_mots,e)
  dist_matrix[ee[1],ee[2]]
})
mintree = mst(g_mots,weights = edgedist)
treedist = sapply(E(mintree),\(e) {
  ee = ends(g_mots,e)
  dist_matrix[ee[1],ee[2]]
})



adj_matrix = dist_matrix<= max(treedist)
g_mots = graph_from_adjacency_matrix(adj_matrix,mode = "undirected",diag = F)
V(g_mots)$name = allnames[1:88]
edge_density(g_mots)
plot(g_mots)
edgedist_ = sapply(E(g_mots),\(e) {
  ee = ends(g_mots,e)
  dist_matrix[ee[1],ee[2]]
})
mg = modularity(g_mots,membership = numcodes,weights = 1/edgedist_)


control_modularity = sapply(1:2000, \(k){
  n = length(numcodes)
  ii = sample(n,n)
  modularity(g_mots,membership = numcodes[ii],weights = 1/edgedist_)
})
hist(control_modularity, main = "Modularity of antibody clusters")
abline(v = mg, col = "red", lwd = 3, lty = 2)
e = ecdf(control_modularity)
e(mg)
# get distance to background----
alldist = compare_motifs(allp_mb,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
colnames(alldist) = allnames
rownames(alldist) = allnames
dist_to_bkg = alldist[,"bkg"]
mean(dist_to_bkg)
which(dist_to_bkg<0.3)
closest = sapply(1:4, \(i){
  v = dist_to_bkg[which(numcodes==i)]
  names(which.min(v))
})

# see which clusters are very similar ----

alldist_sc = scale(alldist)
hist(alldist_sc, breaks = 500, xlim = c(-3,-2))
which(alldist_sc < (-1) & alldist >0, arr.ind = T)
alldist_sc[which(alldist_sc < (-2.4) & alldist >0, arr.ind = T)]

alldist_nobg = alldist[1:91,1:91]
# alldist_nobg[which(alldist_nobg==0)] = NA
alldist_sc = scale(alldist_nobg)
hist(alldist_sc)
hist(alldist_sc, breaks = 500, xlim = c(-3,-2))
which(alldist_sc < (-2.2) & alldist_nobg >0, arr.ind = T)
alldist_sc[which(alldist_sc < (-2.2) & alldist_nobg >0, arr.ind = T)]
