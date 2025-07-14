
# compare kld to kld from the random clusters, calculate z-scores, make graphics
library(stats)
library(gplm)


# get random results
gsizes = sapply(abs,\(ab){
  path = paste0("mixed-7graphs/",ab,"/")
  f = paste0(path,ab,"big7or.RData")
  load(f)
  vcount(G)
})
rnd_results = lapply(gsizes,\(vc) {load(paste0("mixed-7graphs/",vc, "_res_l.RData" ));res})

for(i in 1:4){ 
  data = data.frame(x=log10(rnd_results[[i]][,3]),y=log10(rnd_results[[i]][,1]))
  
  model=kreg(data$x,data$y,grid = data$x)
  ii = order(data$x)
  res = (model$y-data$y[ii])
  plot(model$x,res,col = adjustcolor("black",alpha.f = 0.2),pch = 16,main = paste(abs[i], "log"))
  
}
# calculate z scores from distribution of the kld of random clusters----
kl = all_kld
loglen = log10(lengths)
len = lengths
colcodes = numcodes
z_sc = lapply(1:4,\(i){
  data = data.frame(x=log10(rnd_results[[i]][,3]),y=log10(rnd_results[[i]][,1]))
  real_l = loglen[which(colcodes==i)]
  names(real_l) = names(len[which(colcodes==i)])
  real_k = log10(kl[which(colcodes==i)])
  names(real_k)= names(kl[which(colcodes==i)])
  model=kreg(data$x,data$y,grid = data$x)#Scott's rule of thumb,bandwidth = 11.3
  pred =kreg(data$x,data$y,grid = real_l,bandwidth = model$bandwidth)
  ii = order(data$x)
  sq_res = (model$y-data$y[ii])**2
  #hist(sq_res)
  mod_res = kreg(model$x,sq_res,grid = data$x)#,bandwidth = 11.3
  
  pred_res = kreg(model$x,sq_res,grid = loglen[which(colcodes==i)],bandwidth = mod_res$bandwidth) # var predictions for real data
  sd = sqrt(pred_res$y)
  
  ij = order(real_l)
  
  z = as.vector((real_k[ij]-pred$y)/sd)
  big = which(is.na(pred$y)) # out of the range of the randomly generated clusters
  l = length(model$x)
  
 
  m = which.max(pred$x[-big])
  z[big] = as.vector((real_k[ij][big]-pred$y[-big][m])/sd[-big][m])
  names(z) = names(sort(real_l))
  z
  
})


# p- values
allp = sapply(unlist(z_sc),\(z){
  pnorm(z,lower.tail = F)
})
allp = p.adjust(allp,method = "BH")


# make graphics----
### mapping distances

load(file = "mixed-7graphs/lib_allpeps_7nC.RData")
bkg_m = create_motif(peps_7nC) # motif from library
bgmot = bkg_m@motif


allp_mb = append(allm_wp,bkg_m)  #allp_mb = append(allp_m,bkg_m) for clusters without weights
m_comp_b = compare_motifs(allp_mb,tryRC = FALSE,min.mean.ic = 0.1,method = "EUC")
plot2d = cmdscale(m_comp_b,2)
rownames(plot2d) = c(names(unlist(logpepsets,recursive = F)),"bkg") #rownames(plot2d) = c(names(allpepsets),"bkg")

sign_p  = allp<0.01
# each antibody
nn = nrow(plot2d)
for(i in 1:4){
  plot(plot2d[nn,1],plot2d[nn,2],cex = log((1000))/2,pch = 18,xlim =range( plot2d[,1]),ylim = range(plot2d[,2]),main = paste("clusters of",abs[i]),xlab = "dimension 1",ylab = "dimension 2" )
  
  num = which(colcodes == i)
  o = order(len[num])
  
  p = sign_p[num]
  
  #significant
  points(plot2d[num,][o,][p,],cex = log(sort(len[num][p]))/2,pch = 16,col = adjustcolor("orange", alpha.f = 0.5))
  #non-significant
  points(plot2d[num,][o,],cex = log(sort(len[num]))/2,col = adjustcolor("orange", alpha.f = 0.7))
}

# part of the poster:
# all anti gp120 antibodies together
colors = c("red","green","blue","orange")
plot(plot2d[nn,1],plot2d[nn,2],cex = log((1000))/2,pch = 18,xlim =range( plot2d[,1]),ylim = range(plot2d[,2]),main = paste("Clusters of anti-gp120\nantibodies"),xlab = "dimension 1",ylab = "dimension 2" )

for(i in 2:4){
  
  num = which(colcodes == i)
  o = order(len[num])
  
  p = sign_p[num]
  
  #significant
  points(plot2d[num,][o,][p,],cex = log(sort(len[num][p]))/2,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot2d[num,][o,],cex = log(sort(len[num]))/2,col = adjustcolor(colors[i], alpha.f = 0.7))
}
legend("topright",legend = abs[2:4],col = colors[2:4],pt.cex = 1.6,cex = 1.3,pch = 16,title = "Antibodies",border = F)
legend("bottomleft",x.intersp = 2,pt.cex = log((1000))/2,cex = 1.3,pch = 18,legend = "background")



# only Herceptin
plot(plot2d[nn,1],plot2d[nn,2],cex = log((1000))/2,pch = 18,xlim =range( plot2d[,1]),ylim = range(plot2d[,2]),main = paste("Clusters of Herceptin"),xlab = "dimension 1",ylab = "dimension 2" )

for(i in 1:1){
  
  num = which(colcodes == i)
  o = order(len[num])
  
  p = sign_p[num]
  
  #significant
  points(plot2d[num,][o,][p,],cex = log(sort(len[num][p]))/2,pch = 16,col = adjustcolor(colors[i], alpha.f = 0.5))
  #non-significant
  points(plot2d[num,][o,],cex = log(sort(len[num]))/2,col = adjustcolor(colors[i], alpha.f = 0.7))
}
#legend("topleft",legend = abs[1],col = colors[1],pt.cex = 1.5,pch = 16,title = "Antibodies",border = F)
legend("bottomleft",x.intersp = 2,pt.cex = log((1000))/2,cex = 1.3,pch = 18,legend = "background")
points(plot2d["Herceptin_27",1],plot2d["Herceptin_27",2],pch = 4)


#only legend
plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)

legend("left",legend = c("p>0.01","p<0.01"),col =adjustcolor(colors[1], alpha.f = 0.7),pt.cex = 1.3,pch = c(1,16),title = "p-value",border = F)
legend("center",legend = c("10","50","100","500"),col = adjustcolor(colors[1], alpha.f = 0.7),pch = 16,pt.cex =log(c(10,50,100,500))/2,cex = 1,title = "Number of peptides \nper cluster",x.intersp = 1,y.intersp = 1.2,text.width = 0.15 )

plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
legend("center",legend = abs,col = colors,pt.cex = 1.6,cex = 1.3,pch = 16,title = "Antibodies",border = F)


