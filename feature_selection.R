# select dimensions that best separate the 4 antibodies
require(igraph)
require(uwot)
require(rgl)

load("cmdresult.RData") # codes and cmd_points
load("dist_m.RData") # cluster_dist
codes=c(codes,5)

# Modularity of the graph of clusters by antibody
gGershMabs=graph_from_adjacency_matrix(cluster_dist, mode = "undirected", weighted = T)
ggx=mst(gGershMabs)
# hist(edge_attr(ggx)$weight)
ggx=delete_edges(gGershMabs, E(gGershMabs)[edge_attr(gGershMabs)$weight>0.8])
graph.density(ggx)
ggx=set_edge_attr(ggx, "weight", value=1/edge_attr(ggx)$weight)
modularity(ggx, codes)

# Recursive feature elimination
colnames(cmd_points)=seq_along(cmd_points[1,])
rownames(cmd_points)=seq_along(cmd_points[,1])
x=rfe(t(cmd_points),codes[-length(codes)])
bestff=rownames(x[[1]])

X=cmd_points[,bestff]
ux=umap(X, n_neighbors = 20, min_dist = 0.1, spread=7,  n_epochs = 1500,
        verbose=T, n_components = 3, init = "normlaplacian",n_sgd_threads=14)
oldpal=palette()
palette(c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),rgb(0,0,1,0.5),"#FFA50090"))
plot(ux, col=codes[-length(codes)], cex=2, pch=16)
legend("bottomleft", legend=c("Herceptin","21c","17b","b12"), fill=1:4)
colnames(ux)=c("D1","D2","D3")
plot3d(ux, col=codes[-length(codes)], size=10)
palette(oldpal)
