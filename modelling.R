args<-commandArgs(TRUE);

data<-scan(file=paste('data_sample/sample',args[1],args[2],'txt',sep='.'),sep='\t',what=numeric(0));
data<-matrix(data,2535,1000*as.numeric(as.character(args[1])));
#data<-matrix(data,2535,length(data)/2535);
data<-t(data);

#data<-scan(file='data_sample/sample.1000.1.notab.txt',what=character(0));
#data<-strsplit(data,split='');
#data<-t(as.data.frame(data));
#mode(data)<-"numeric";
#rownames(data)<-NULL;

#
args<-c(1000,1);
inds_all<-read.table(file='merged_samples_lj.20180712.ALL.panel',header=TRUE);
inds_all[,1]<-as.character(inds_all[,1]);
inds_all[,2]<-as.character(inds_all[,2]);
inds_all[,3]<-as.character(inds_all[,3]);
inds_all[,4]<-as.character(inds_all[,4]);

inds_super_pop_all<-list();
inds_pop_all<-list();
inds_super_pop_names_all<-unique(inds_all[,3]);
inds_pop_names_all<-unique(inds_all[,2]);

super_pop2pop<-unique(inds_all[,2:3]);
super_pop_map<-list();
pop_map<-numeric(length(inds_pop_names_all));


for(i in 1:length(inds_super_pop_names_all)){
	inds_super_pop_all[[i]]<-which(inds_all[,3]==inds_super_pop_names_all[i]);
	super_pop_map[[i]]<-match(super_pop2pop[super_pop2pop[,2]==inds_super_pop_names_all[i],1],inds_pop_names_all);
}
for(i in 1:length(inds_pop_names_all)){
	inds_pop_all[[i]]<-which(inds_all[,2]==inds_pop_names_all[i]);
	pop_map[i]<-which(inds_super_pop_names_all==super_pop2pop[super_pop2pop[,1]==inds_pop_names_all[i],2]);
}

less_sample<-1;
inds_pop_less<-list();
inds_super_pop_less<-list();
for(i in 1:length(inds_pop_names_all)){
	inds_pop_less[[i]]<-sample(inds_pop_all[[i]],round(length(inds_pop_all[[i]])*less_sample));
}
for(i in 1:length(inds_super_pop_names_all)){
	inds_super_pop_less[[i]]<-do.call(c,inds_pop_less[super_pop_map[[i]]]);
}
ind_idx_less<-do.call(c,inds_pop_less);
ind_idx_less<-ind_idx_less[order(ind_idx_less)];
ind_num_less<-length(ind_idx_less);

inds_super_pop<-inds_super_pop_less;
inds_pop<-inds_pop_less;
inds_super_pop_names<-inds_super_pop_names_all;
inds_pop_names<-inds_pop_names_all;

inds_less<-inds_all[ind_idx_less,];
data<-data[,ind_idx_less];

rownames(data)<-NULL;
colnames(data)<-NULL;
var_num<-dim(data)[1];
ind_num<-dim(data)[2];

ind_idx<-1:ind_num;
var_idx<-1:var_num;

library(slam);
library(topicmodels);

find_word <- function(doc){
	which(doc!=0);
}
js<-apply(data,2,find_word);
j<-do.call(c,js);
ns<-sapply(js,length);
i<-rep(1:ind_num,ns);
v<-data[data!=0];
x<-simple_triplet_matrix(i,j,v,nrow=ind_num,ncol=var_num,dimnames=NULL);

times<-list();

times[[1]]<-Sys.time();
times[[1]];
lda<-LDA(x,k=17,method="Gibbs",control=list(iter=500));
times[[2]]<-Sys.time();
times[[2]];

results<-posterior(lda);

delta<-distHellinger(results$topics);

#rm(data,x,i,j,js,ns,v);
save(inds_super_pop,inds_pop,inds_super_pop_names,inds_pop_names,times,lda,results,delta,ind_idx_less,inds_less,ind_num_less,file=paste('sample_models/sample',args[1],args[2],'RData',sep='.'));

q(save='no');
