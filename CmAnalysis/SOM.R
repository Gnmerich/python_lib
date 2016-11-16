#!/bin/Rscript

library("optparse")
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, help="CM feature vectors", metavar="character"),
    make_option(c("-x", "--xdim"), type="integer", default=5, help="Dimension X of SOM grid", metavar="integer" ),
    make_option(c("-y", "--ydim"), type="integer", default=5, help="Dimension Y of SOM grid", metavar="integer" ),
    make_option(c("-r", "--rlen"), type="integer", default=500, help="SOM Learning iterations", metavar="integer" ),
	make_option(c("-c", "--clu"), type="integer", default=15, help="desired number of clusters", metavar="integer"),
	make_option(c("-o", "--out"), type="character", default="out.txt", help="output file prefix", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("No input file supplied!", call.=FALSE)
}

####################################################################################################################

require(kohonen)
require(graphics)

#Create SOM
data<-read.table(file=opt$file, header=TRUE, row.names=1)
data_train_matrix<-as.matrix(data)
som_grid <- somgrid(xdim=opt$xdim, ydim=opt$ydim, topo="rectangular")
som_model <- som(data_train_matrix, grid=som_grid, rlen=opt$rlen, alpha=c(0.05,0.01), keep.data=TRUE, n.hood='circular')

#Clustering - find 'optimal' number with kmeans
mydata <- som_model$codes
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:(opt$xdim*opt$ydim-1)) {
    wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}

#Plot Kmeans results TODO: use actual kmeans result in clustering
pdf(file='SOM_kmeans.pdf')
plot(wss)
dev.off()

#Hierarchical clustering
som_cluster <- cutree(hclust(dist(som_model$codes)), opt$clu)
pretty_palette<-palette(rainbow(opt$clu))
dev.off()
pretty_palette<-palette(rainbow(opt$clu))
#dev.off()

#Plot clustered SOM
pdf(file='SOM_clustered.pdf')
plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters")
add.cluster.boundaries(som_model, som_cluster)
dev.off()

#Plot SOM training process
pdf(file='SOM_changes.pdf')
plot(som_model, type="changes")
dev.off()
#Plot SOM with feature vectors
pdf(file='SOM.pdf')
plot(som_model, type='codes')
dev.off()
#Plot SOM with feature vectors and cluster boundaries
pdf(file='SOM_bound.pdf')
plot(som_model, type='codes')
add.cluster.boundaries(som_model, som_cluster)
dev.off()
#Plot SOM with pop densities
pdf(file='SOM_count.pdf')
plot(som_model, type="count")
dev.off()
#Write each bin's vector to file
write.table(som_model$codes, file='SOM_codes', sep='\t')
write.table(som_model$data, file='SOM_data', sep='\t')


#Write Mapping Vector->bin vector->cluster to file
vec_names<-rownames(som_model$data)
sink(file='SOM_clusters')
cat(paste('CLUSTER','\t','BIN','\t','NAME','\t'))
cat(paste(sep='\t', colnames(som_model$data)),'\n')
for (i in 1:length(som_model$unit.classif)) {
    cat(som_cluster[som_model$unit.classif[i]])
    cat('\t')
    cat(som_model$unit.classif[i])
    cat('\t')
    cat(vec_names[i])
    cat('\t')
    cat(som_model$data[i,])
    cat('\n')
}
sink()
