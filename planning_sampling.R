#initialize var list
var<-list();
for(i in 1:23){
	if(i<23){
		var[[i]]<-read.table(file=paste('1000g_slim/variants.chr',i,'.txt',sep=''),sep='\t',header=FALSE);
	}else{
		var[[i]]<-read.table(file='1000g_slim/variants.chrX.txt',sep='\t',header=FALSE);
	}
	var[[i]][,1]<-as.character(var[[i]][,1]);
	var[[i]][,2]<-as.numeric(as.character(var[[i]][,2]));
	var[[i]][,3]<-as.character(var[[i]][,3]);
	var[[i]][,4]<-as.character(var[[i]][,4]);
	var[[i]][,5]<-as.character(var[[i]][,5]);
	var[[i]][,6]<-as.numeric(as.character(var[[i]][,6]));
	var[[i]][,7]<-as.numeric(as.character(var[[i]][,7]));
	var[[i]][,8]<-as.numeric(as.character(var[[i]][,8]));
	var[[i]][,9]<-as.numeric(as.character(var[[i]][,9]));
	var[[i]][,10]<-as.numeric(as.character(var[[i]][,10]));
	var[[i]][,11]<-as.numeric(as.character(var[[i]][,11]));
	var[[i]][,12]<-as.numeric(as.character(var[[i]][,12]));
}

chr_vars<-sapply(var,dim)[1,];
af_boundaries<-cbind(rep(0,6),rep(1,6));

af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.1;

filtered_idx<-list();
for(i in 1:23){
	filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
							var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
							var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
							var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
							var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
							var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
}
filtered_idx<-do.call(rbind,filtered_idx);

plan_schema<-c(1000,2000,3000,5000,7000,10000,15000,20000,30000,50000,70000,100000,150000,200000,300000,500000,1000000);
rep_num<-c(rep(10,15),1,1);
plans<-matrix(0,sum(plan_schema*rep_num),4);
for(i in 1:length(plan_schema)){
	for(j in 1:rep_num[i]){
		idx<-sample(1:dim(filtered_idx)[1],plan_schema[i]);
		range<-(sum(plan_schema[0:(i-1)]*rep_num[0:(i-1)])+plan_schema[i]*(j-1)+1):(sum(plan_schema[0:(i-1)]*rep_num[0:(i-1)])+plan_schema[i]*j);
		plans[range,]<-cbind(plan_schema[i]/1000,j,filtered_idx[idx,]);
	}
}
plans<-plans[order(plans[,3],plans[,4],plans[,1],plans[,2]),];
write.table(plans,file='sampling_plans_01.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

#find best af region
plan_schema<-500000;
afs<-c(0.1,0.125,0.15,0.175,0.2,0.225,0.25);
rep_num<-20;

plans_af3<-matrix(0,plan_schema*length(afs)*rep_num,4);
filtered_idx<-list();
plans_af3_map<-matrix(0,length(afs)*rep_num,4);

for(afi in 1:length(afs)){
	af_boundaries[1,1]<-0.001;
	af_boundaries[1,2]<-afs[afi];
	filtered_idx<-list();
	for(i in 1:23){
		filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
										var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
										var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
										var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
										var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
										var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
	}
	filtered_idx<-do.call(rbind,filtered_idx);
	for(j in 1:rep_num){
		idx<-sample(1:dim(filtered_idx)[1],plan_schema);
		range<-((afi-1)*rep_num*plan_schema+(j-1)*plan_schema+1):((afi-1)*rep_num*plan_schema+j*plan_schema);
		plans_af3[range,]<-cbind(plan_schema/1000,afi*rep_num+j,filtered_idx[idx,]);
		plans_af3_map[(afi-1)*rep_num+j,]<-c(plan_schema/1000,afi*rep_num+j,0.001,afs[afi]);
	}
}
plans_af3<-plans_af3[order(plans_af3[,3],plans_af3[,4],plans_af3[,1],plans_af3[,2]),];
write.table(plans_af3,file='sampling_plans_af3.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
write.table(plans_af3_map,file='sampling_plans_af3_map.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.5;

#different af
plan_schema<-100000;
afs<-seq(0,0.45,0.05);
rep_num<-10;

plans_af<-matrix(0,plan_schema*length(afs)*rep_num,4);
filtered_idx<-list();
plans_af_map<-matrix(0,length(afs)*rep_num,4);
for(afi in 1:length(afs)){
	af_boundaries[1,1]<-max(afs[afi],0.001);
	af_boundaries[1,2]<-afs[afi]+0.05;
	filtered_idx<-list();
	for(i in 1:23){
		filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
										var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
										var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
										var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
										var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
										var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
	}
	filtered_idx<-do.call(rbind,filtered_idx);
	for(j in 1:rep_num){
		idx<-sample(1:dim(filtered_idx)[1],plan_schema);
		range<-((afi-1)*rep_num*plan_schema+(j-1)*plan_schema+1):((afi-1)*rep_num*plan_schema+j*plan_schema);
		plans_af[range,]<-cbind(plan_schema/1000,afi*rep_num+j,filtered_idx[idx,]);
		plans_af_map[(afi-1)*rep_num+j,]<-c(plan_schema/1000,afi*rep_num+j,max(afs[afi],0.001),afs[afi]+0.05);
	}
}
plans_af<-plans_af[order(plans_af[,3],plans_af[,4],plans_af[,1],plans_af[,2]),];
write.table(plans_af,file='sampling_plans_af.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
write.table(plans_af_map,file='sampling_plans_af_map.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.5;

#different regions

#initialize genome information
chr_length<-read.table(file='~/genomes/hg19/hg19.fa.fai',sep='\t',header=FALSE);
chr_length<-chr_length[1:23,2];
options(scipen=20);
span<-10000000;
genome_divisions<-list();
g_idx<-1;
for(i in 1:23){
	for(regions in seq(1,chr_length[i]-span/2,span)){
		genome_divisions[[g_idx]]<-c(i,regions,regions+span-1,length(which(var[[i]][,2]<regions))+1,length(which(var[[i]][,2]<regions+span-1)));
		g_idx<-g_idx+1;
	}
	genome_divisions[[g_idx-1]][3]<-chr_length[i];
	genome_divisions[[g_idx-1]][5]<-chr_vars[i];
}
genome_divisions<-t(as.data.frame(genome_divisions));
rownames(genome_divisions)<-NULL;

#select variants
af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.5;

filtered_idx<-list();
for(i in 1:23){
	filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
							var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
							var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
							var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
							var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
							var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
}
division_idx<-list();
for(i in 1:dim(genome_divisions)[1]){
	idx<-which(filtered_idx[[genome_divisions[i,1]]][,2]>=genome_divisions[i,4]&filtered_idx[[genome_divisions[i,1]]][,2]<=genome_divisions[i,5]);
	if(length(idx)>100){
		division_idx[[i]]<-cbind(0,i,filtered_idx[[genome_divisions[i,1]]][idx,]);
	} else {
		division_idx[[i]]<-matrix(0,0,4);
	}
}

plans_division<-do.call(rbind,division_idx);
plans_division<-plans_division[order(plans_division[,3],plans_division[,4],plans_division[,1],plans_division[,2]),];
plans_division_map<-cbind(0,1:dim(genome_divisions)[1],genome_divisions,sapply(division_idx,dim)[1,]);
write.table(plans_division,file='sampling_plans_divisions.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
write.table(plans_division_map,file='sampling_plans_divisions_map.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

#different variant density

af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.5;
plan_schema<-100000;
rep_num<-10;
block_schema<-c(1,2,3,5,7,10,15,20,30,50,70,100,150,200,250);
plans_density<-matrix(0,length(block_schema)*rep_num*plan_schema,4);
plans_density_map<-matrix(0,length(block_schema)*rep_num,7);
filtered_idx<-list();
for(i in 1:23){
	filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
							var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
							var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
							var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
							var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
							var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
}
filtered_idx<-do.call(rbind,filtered_idx);
for(block in 1:length(block_schema)){
	blocks<-sample(1:(dim(filtered_idx)[1]-plan_schema*block_schema[block]),rep_num);
	for(j in 1:rep_num){
		range<-((block-1)*rep_num*plan_schema+(j-1)*plan_schema+1):((block-1)*rep_num*plan_schema+j*plan_schema);
		idx<-seq(blocks[j],(blocks[j]+plan_schema*block_schema[block]-1),block_schema[block]);
		plans_density[range,]<-cbind(plan_schema/1000,block*rep_num+j+100,filtered_idx[idx,]);

		chr_switch<-which(filtered_idx[idx[1:(plan_schema-1)],1]-filtered_idx[idx[2:plan_schema],1]!=0); #chr span idx in sampled filtered_idx
		chr_switch<-cbind(c(1,chr_switch+1),c(chr_switch,plan_schema)); #start idx and end idx for each chr in sampled filtered_idx
		avg_distance<-0;
		for(k in 1:dim(chr_switch)[1]){
			start_chr<-filtered_idx[idx[chr_switch[k,1]],1];#start_chr == end_chr, guarenteed by previous step
			start_idx<-filtered_idx[idx[chr_switch[k,1]],2];
			end_idx<-filtered_idx[idx[chr_switch[k,2]],2];
			avg_distance<-avg_distance+var[[start_chr]][end_idx,2]-var[[start_chr]][start_idx,2];
		}
		avg_distance<-avg_distance/(plan_schema-dim(chr_switch)[1]);
		plans_density_map[(block-1)*rep_num+j,]<-c(plan_schema/1000,block*rep_num+j+100,block_schema[block],j,blocks[j],blocks[j]+plan_schema-1,avg_distance);
	}
}
plans_density<-plans_density[order(plans_density[,3],plans_density[,4],plans_density[,1],plans_density[,2]),];
write.table(plans_density,file='sampling_plans_density.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
write.table(plans_density_map,file='sampling_plans_density_map.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

#different variant type (snv/ins/del)
af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.5;

plan_schema<-100000;
rep_num<-10;

filtered_snv_idx<-list();
filtered_indel_idx<-list();

plans_type<-matrix(0,2*rep_num*plan_schema,4);
for(i in 1:23){
	filtered_snv_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
							var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
							var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
							var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
							var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
							var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]&
							nchar(var[[i]][,4])==nchar(var[[i]][,5])));
	filtered_indel_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
							var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
							var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
							var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
							var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
							var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]&
							nchar(var[[i]][,4])!=nchar(var[[i]][,5])));
}
filtered_snv_idx<-do.call(rbind,filtered_snv_idx);
filtered_indel_idx<-do.call(rbind,filtered_indel_idx);

for(j in 1:rep_num){
	range<-((j-1)*plan_schema+1):(j*plan_schema);
	idx<-sample(1:dim(filtered_snv_idx)[1],plan_schema);
	plans_type[range,]<-cbind(plan_schema/1000,j+260,filtered_snv_idx[idx,]);
}
for(j in 1:rep_num){
	range<-plan_schema*rep_num+((j-1)*plan_schema+1):(j*plan_schema);
	idx<-sample(1:dim(filtered_indel_idx)[1],plan_schema);
	plans_type[range,]<-cbind(plan_schema/1000,j+270,filtered_indel_idx[idx,]);
}
plans_type<-plans_type[order(plans_type[,3],plans_type[,4],plans_type[,1],plans_type[,2]),];
write.table(plans_type,file='sampling_plans_type.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

#####for stat these variants and plot
snp_idx<-unique(plans_type[plans_type[,2]<=270,3:4]);
indel_idx<-unique(plans_type[plans_type[,2]>270,3:4])
snps<-var[[1]][snp_idx[snp_idx[,1]==1,2],];
indels<-var[[1]][indel_idx[indel_idx[,1]==1,2],];
for(i in 2:23){
	snps<-rbind(snps,var[[i]][snp_idx[snp_idx[,1]==i,2],]);
	indels<-rbind(indels,var[[i]][indel_idx[indel_idx[,1]==i,2],]);
}
c<-list();
for(i in 1:6){
	c[[(i-1)*2+1]]<-snps[snps[,i+5]!=0,i+5];
	c[[i*2]]<-indels[indels[,i+5]!=0,i+5];
}
bars<-c("SNP_ALL","INDEL_ALL","SNP_EAS","INDEL_EAS","SNP_AMR","INDEL_AMR","SNP_AFR","INDEL_AFR","SNP_EUR","SNP_INDEL","SNP_SAS","INDEL_SAS");
d<-as.data.frame(matrix(0,0,2));
for(i in 1:6){
	d<-rbind(d,cbind(bars[(i-1)*2+1],snps[snps[,i+5]!=0,i+5]));
	d<-rbind(d,cbind(bars[i*2],indels[indels[,i+5]!=0,i+5]));
}
d[,1]<-as.factor(d[,1]);

redc<-rgb(251,128,114, maxColorValue = 255);
bluec<-rgb(128,177,211, maxColorValue = 255);
greyc<-rgb(89,89,89, maxColorValue = 255);
png('snp_indel_boxplot.png',2000,1200);
boxplot(c,varwidth=TRUE,ylim=c(0,0.5),col=rep(c(redc,bluec),6),border=rep(greyc,12),names=bars,lwd=3,pch=20);
dev.off();

#use coding variants only
af_boundaries[1,1]<-0.001;
af_boundaries[1,2]<-0.5;

filtered_idx<-list();
for(i in 1:23){
	filtered_idx[[i]]<-cbind(i,which(var[[i]][,12]==1&
							var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
							var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
							var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
							var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
							var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
							var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
}
filtered_idx<-do.call(rbind,filtered_idx);
plans_coding<-cbind(round(dim(filtered_idx)[1]/1000),1,filtered_idx);
plans_coding<-plans_coding[order(plans_coding[,3],plans_coding[,4],plans_coding[,1],plans_coding[,2]),];
write.table(plans_coding,file='sampling_plans_coding.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

##different af in each super pop

plans_af_and_pop<-matrix(0,600000,4);
plan_idx<-1;
af_boundaries[1,1]<-0;
af_boundaries[1,2]<-1;
for(af_idx in 2:6){
	for(af_value in seq(0,0.25,0.05)){
		af_boundaries[af_idx,1]<-max(af_value,0.001);
		af_boundaries[af_idx,2]<-af_value+0.05;
		filtered_idx<-list();
		for(i in 1:23){
			filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
											var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
											var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
											var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
											var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
											var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
		}
		filtered_idx<-do.call(rbind,filtered_idx);
		idx<-sample(1:dim(filtered_idx)[1],20000);
		plans_af_and_pop[((plan_idx-1)*20000+1):(plan_idx*20000),]<-cbind(20,plan_idx+15,filtered_idx[idx,]);
		plan_idx<-plan_idx+1;
	}
	af_boundaries[af_idx,1]<-0;
	af_boundaries[af_idx,2]<-1;
}
plans_af_and_pop<-plans_af_and_pop[order(plans_af_and_pop[,3],plans_af_and_pop[,4],plans_af_and_pop[,1],plans_af_and_pop[,2]),];
write.table(plans_af_and_pop,file='sampling_plans_af_and_pop.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
plan_map<-cbind(16:(plan_idx+14),super_pop[gl(5,6)],rep(seq(0,0.25,0.05),5),rep(seq(0,0.25,0.05),5)+0.05);
write.table(plan_map,file='sampling_plans_af_and_pop_map.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

#######
#different variant density x different af
afs<-seq(0,0.45,0.05);
plan_schema<-100000;
af_varnums<-c(21833853,1730306,1061003,790497,647611,543811,468352,414062,364646,330163);
block_schema_percent<-c(0.01,0.03,0.06,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1);
block_schema_start<-c(1,3,4,5,6,6,7,7,8,9);
plans_density<-matrix(0,sum(length(block_schema_percent)-block_schema_start+1)*plan_schema,4);
plans_density_map<-matrix(0,sum(length(block_schema_percent)-block_schema_start+1),8);

range_start<-1;
for(afi in 1:length(afs)){
	af_boundaries[1,1]<-max(afs[afi],0.001);
	af_boundaries[1,2]<-afs[afi]+0.05;
	filtered_idx<-list();
	for(i in 1:23){
		filtered_idx[[i]]<-cbind(i,which(var[[i]][,6]>=af_boundaries[1,1]&var[[i]][,6]<=af_boundaries[1,2]&
								var[[i]][,7]>=af_boundaries[2,1]&var[[i]][,7]<=af_boundaries[2,2]&
								var[[i]][,8]>=af_boundaries[3,1]&var[[i]][,8]<=af_boundaries[3,2]&
								var[[i]][,9]>=af_boundaries[4,1]&var[[i]][,9]<=af_boundaries[4,2]&
								var[[i]][,10]>=af_boundaries[5,1]&var[[i]][,10]<=af_boundaries[5,2]&
								var[[i]][,11]>=af_boundaries[6,1]&var[[i]][,11]<=af_boundaries[6,2]));
	}
	filtered_idx<-do.call(rbind,filtered_idx);
	block_schema<-floor(dim(filtered_idx)[1]/plan_schema*block_schema_percent[block_schema_start[afi]:length(block_schema_percent)]*1000)/1000;

	for(block in 1:length(block_schema)){
		blocks<-sample(1:(dim(filtered_idx)[1]-plan_schema*block_schema[block]),1);
		range<-range_start:(range_start+plan_schema-1);
		range_start<-range_start+plan_schema;
		experiment_idx<-sum(length(block_schema_percent)-block_schema_start[0:(afi-1)]+1)+block;

		idx<-floor(seq(blocks,(blocks+plan_schema*block_schema[block]-1),block_schema[block]));
		plans_density[range,]<-cbind(plan_schema/1000,800+experiment_idx,filtered_idx[idx,]);

		chr_switch<-which(filtered_idx[idx[1:(plan_schema-1)],1]-filtered_idx[idx[2:plan_schema],1]!=0); #chr span idx in sampled filtered_idx
		chr_switch<-cbind(c(1,chr_switch+1),c(chr_switch,plan_schema)); #start idx and end idx for each chr in sampled filtered_idx
		avg_distance<-0;
		for(k in 1:dim(chr_switch)[1]){
			start_chr<-filtered_idx[idx[chr_switch[k,1]],1];#start_chr == end_chr, guarenteed by previous step
			start_idx<-filtered_idx[idx[chr_switch[k,1]],2];
			end_idx<-filtered_idx[idx[chr_switch[k,2]],2];
			avg_distance<-avg_distance+var[[start_chr]][end_idx,2]-var[[start_chr]][start_idx,2];
		}
		avg_distance<-avg_distance/(plan_schema-dim(chr_switch)[1]);
		plans_density_map[experiment_idx,]<-c(plan_schema/1000,800+experiment_idx,block_schema[block],blocks,blocks+plan_schema*block_schema[block]-1,avg_distance,af_boundaries[1,1],af_boundaries[1,2]);
	}
}
plans_density<-plans_density[order(plans_density[,3],plans_density[,4],plans_density[,1],plans_density[,2]),];
write.table(plans_density,file='sampling_plans_densityXaf.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
write.table(plans_density_map,file='sampling_plans_densityXaf_map.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);

##save and quit

super_pop<-c("EAS","AMR","AFR","EUR","SAS");
save.image(file='var.RData');
q(save='no');
