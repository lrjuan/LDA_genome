args<-commandArgs(TRUE);

data<-scan(file=paste('data_sample/sample',args[1],args[2],'txt',sep='.'),sep='\t',what=numeric(0));
data<-matrix(data,2535,1000*as.numeric(as.character(args[1])));
data<-t(data);

hamming<-matrix(0,2535,2535);
weight<-numeric(2535);

for(i in 1:dim(data)[1]){
	idx<-which(data[i,]==0);
	hamming[idx,-idx]<-hamming[idx,-idx]+1;
	hamming[-idx,idx]<-hamming[-idx,idx]+1;
	weight[-idx]<-weight[-idx]+1;
}

save(hamming,weight,file=paste('sample_models/sample.hamming',args[1],args[2],'RData',sep='.'));

q(save='no');
