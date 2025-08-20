require(FactoMineR)
require(factoextra)
require(protr)

pcaaind=PCA(scale(t(AAindex[,7:26])), ncp=19)
fviz_screeplot(pcaaind, ncp=19)
aaindnew=pcaaind$ind$coord[,1:11]
maaind=as.matrix(dist(aaindnew))

distm_ = cbind(maaind,0)
colnames(distm_)[21] = "-"

# convert to similarity matrix
simm_ = max(distm_) - distm_
simm_ = (simm_ - min(simm_))/diff(range(simm_))

#try with blossum weights
data("BLOSUM45")
blosum = BLOSUM45[1:20,]
blosum = cbind(blosum,rep(0,20))
colnames(blosum)[ncol(blosum)] = "-"

# contact correlations from paper 

contm = read.table("aa_contacts.txt",header = T)[,1:20]
contm = as.matrix(contm)
contm_ = scale(contm)
corrmat = cor(t(contm_))
heatmap(corrmat,scale = "none")

# zScales scores

library(Peptides)
data("AAdata")
zmatrix = do.call(cbind,AAdata$zScales)
zcorr = cor(t(zmatrix))
heatmap(zcorr,scale = "none")
