
load(file = "epitope_graphs.RData")
load(file = "mixed-7graphs/ppm_pssm_w.RData")
for(i in 1:4){
  names(all_ppm[[i]]) = clnames[[i]]
}
ppms = unlist(all_ppm,recursive = F)

require(ggseqlogo)
require(ggplot2)
require(Biostrings)


library(vctrs)
#try with blossum weights
data("BLOSUM45")
blosum = BLOSUM45[1:20,]
blosum = cbind(blosum,rep(0,20))
colnames(blosum)[ncol(blosum)] = "-"

max_scores_all = future_sapply(1:4,\(ab_i){
  paths = longpaths(g_full[[ab_i]])
  sapply(ppms,\(mot){
    
    scores = unlist(sapply(paths,\(paths2) {
      paths2 = list_drop_empty(paths2)
      sapply(paths2,\(path){
        n = length(path)
        if(n<=ncol(mot)){
          max(sapply(1:(ncol(mot)-n+1),\(i){
            #sum(diag(mot[epi,i:(i+n-1)]))
            m = mot[,i:(i+n-1)]*cbind(blosum[,path]) # part of motif times substitution frequency of path aas
            mean(colSums(m))
          }))
        } 
        else{
          max(sapply(1:(n - ncol(mot)+1),\(i){
            #sum(diag(mot[epi[i:(i+ncol(mot)-1)],]))
            m = mot*cbind(blosum[,path[i:(i+ncol(mot)-1)]])
            mean(colSums(m))
          }))
        }
      })
    }))
    max(scores)
  })
})

# plot score distributions
library(vioplot)
colors = c("red","green","blue","orange")
for(i in 1:4){
  l = sapply(logfreqs[[i]],sum)
  o = order(l)
  sc = max_scores_all[which(codes==abs[i]),i][o]
  control = max_scores_all[,i]
  e = ecdf(control)
  qt = e(sc)
  qth = which(qt>0.9)
  plot(sc, col = colors[i], ylim = range(max_scores_all),pch = 16,cex = 1.3,main = abs[i],ylab = "scores",xlab = "smallest to largest cluster")
  
  colin = "lightgrey"
  colbd = "grey"
  vioplot(control,add = T,at = 1.5,wex = 1.5, col=colin, border=colbd, rectCol=colbd, lineCol=colbd)
  points(sc, col = colors[i],pch = 16,cex = 1.3)
  points(sc, cex = 1.3)
  if(max(qt)>0.9){
    x = (1:length(sc)+0.2)
    y = (sc+ sd(max_scores_all)/4)
    text(x[qth], y[qth], labels = round(qt[qth],2), cex = 1.2)
    
  }
  
}
