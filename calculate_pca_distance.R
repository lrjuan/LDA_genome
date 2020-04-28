args<-commandArgs(TRUE);

data<-scan(file=paste('data_sample/sample',args[1],args[2],'txt',sep='.'),sep='\t',what=numeric(0));
data<-matrix(data,2535,length(data)/2535);
#data<-t(data);

#pc<-princomp(data,cor=TRUE);
times<-list();

times[[1]]<-Sys.time();
times[[1]];
pc<-prcomp(data);
times[[2]]<-Sys.time();
times[[2]];
dist_matrice<-list();
for(i in 1:100){
		dist_matrice[[i]]<-as.matrix(dist(pc$x[,1:i]));
}
point_coordinates<-pc$x[,1:2];
dist8_matrix<-as.matrix(dist(pc$x[,1:8]));
dist17_matrix<-as.matrix(dist(pc$x[,1:17]));
dist50_matrix<-as.matrix(dist(pc$x[,1:50]));
dist_all_matrix<-as.matrix(dist(pc$x));
pcrotation<-pc$rotation[,1:17];
pcx17<-pc$x[,1:17];
impor17<-summary(pc)$importance[,1:17];
save(pc,file=paste('sample_models/sample.pca',args[1],args[2],'RData',sep='.'));
save(dist_matrice,file=paste('sample_models/sample.pca_dist1_100',args[1],args[2],'RData',sep='.'));
save(pcrotation,dist_all_matrix,dist8_matrix,dist17_matrix,dist50_matrix,point_coordinates,pcx17,impor17,file=paste('sample_models/sample.pca_slim',args[1],args[2],'RData',sep='.'));

q(save='no');
