#lda-mds

load('lda.RData')

male_list<-which(inds[,4]=="male");
female_list<-which(inds[,4]=="female");

plot_scale=1200;
png('all_conf.png',plot_scale,plot_scale);
#colors<-c("#1F78B4","#E31A1C","#33A02C","#FDBF6F","#FB9A99","#A6CEE3","#B2DF8A");
colors<-c("#377EB8","#E41A1C","#FF7F00","#984EA3","#4DAF4A","#A65628","#FFFF33");
pchs<-21:25;
plot(all_conf[,1],all_conf[,2],type='n',pch=16);
for(i in 1:5){
	points(all_conf[inds_super_pop[[i]],1],all_conf[inds_super_pop[[i]],2],pch=16,col=colors[i],cex=2);
}
dev.off();

png('all_conf_pop.png',plot_scale,plot_scale);
plot(all_conf[,1],all_conf[,2],type='n',pch=16);
for(i in 1:5){
	for(j in 1:length(inds_pop2super[[i]])){
		k<-inds_pop[[inds_pop2super[[i]][j]]];
		points(all_conf[k,1],all_conf[k,2],pch=(20+i),col=colors[j],bg=colors[j],cex=2);
	}
}
dev.off();

for(i in 1:5){
	png(paste(inds_super_pop_names[[i]],'pop.png',sep='_'),plot_scale,plot_scale);
	plot(inds_super_pop_conf[[i]][,1],inds_super_pop_conf[[i]][,2],type='n',pch=16);
	for(j in 1:length(inds_pop2super[[i]])){
		k<-match(inds_pop[[inds_pop2super[[i]][j]]],inds_super_pop[[i]]);
		points(inds_super_pop_conf[[i]][k,1],inds_super_pop_conf[[i]][k,2],pch=16,col=colors[j],cex=3);
#		points(inds_super_pop_conf[[i]][k,1],inds_super_pop_conf[[i]][k,2],pch=(20+i),col=colors[j],bg=colors[j])
	}
	dev.off();
#	png(paste(inds_super_pop_names[[i]],'gender.png',sep='_'),plot_scale,plot_scale);
#	plot(inds_super_pop_conf[[i]][,1],inds_super_pop_conf[[i]][,2],type='p',pch=20);
#	k<-match(male_list,inds_super_pop[[i]]);
#	k<-k[!is.na(k)];
#	points(inds_super_pop_conf[[i]][k,1],inds_super_pop_conf[[i]][k,2],pch=20,col=colors[1]);
#	k<-match(female_list,inds_super_pop[[i]]);
#	k<-k[!is.na(k)];
#	points(inds_super_pop_conf[[i]][k,1],inds_super_pop_conf[[i]][k,2],pch=20,col=colors[2]);
#	dev.off();
}

#for(i in 1:26){
#	png(paste(inds_pop_names[[i]],'gender.png',sep='_'),plot_scale,plot_scale);
#	plot(inds_pop_conf[[i]][,1],inds_pop_conf[[i]][,2],type='p',pch=20);
#	k<-match(male_list,inds_pop[[i]]);
#	k<-k[!is.na(k)];
#	points(inds_pop_conf[[i]][k,1],inds_pop_conf[[i]][k,2],pch=20,col=colors[1]);
#	k<-match(female_list,inds_pop[[i]]);
#	k<-k[!is.na(k)];
#	points(inds_pop_conf[[i]][k,1],inds_pop_conf[[i]][k,2],pch=20,col=colors[2]);
#	dev.off();
#}

q(save='no');
