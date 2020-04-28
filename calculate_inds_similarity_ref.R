pop_relations<-read.table(file='pop_relation.txt',sep='\t',header=TRUE);
pop_names<-names(pop_relations);
pop_relations<-as.matrix(pop_relations);
rownames(pop_relations)<-colnames(pop_relations);


#inds<-read.table(file='integrated_call_samples_v3.20130502.ALL.panel',header=TRUE);
inds<-read.table(file='merged_samples_lj.20180712.ALL.panel',header=TRUE);
inds[,1]<-as.character(inds[,1]);
inds[,2]<-as.character(inds[,2]);
inds[,3]<-as.character(inds[,3]);
inds[,4]<-as.character(inds[,4]);
inds<-cbind(inds,0);
inds[,5]<-match(inds[,2],pop_names);

ind_num<-dim(inds)[1];

inds_similarity<-matrix(9,ind_num,ind_num);

t<-cbind(rep(1:ind_num,rep(ind_num,ind_num)),rep(1:ind_num,ind_num));
inds_similarity[t]<-pop_relations[cbind(inds[t[,1],5],inds[t[,2],5])];
rownames(inds_similarity)<-inds[,1];
colnames(inds_similarity)<-inds[,1];

##relatives included
inds_relations<-read.table(file='inds_relations.txt',sep='\t');
inds_relations[,1]<-as.character(inds_relations[,1]);
inds_relations[,2]<-as.character(inds_relations[,2]);
inds_relations[,3]<-as.character(inds_relations[,3]);
idx<-match(inds_relations[,1],inds[,1]);
inds_relations<-inds_relations[!is.na(idx),];
idx<-match(inds_relations[,2],inds[,1]);
inds_relations<-inds_relations[!is.na(idx),];
idx<-which(inds_relations[,1]==inds_relations[,2]);
inds_relations<-inds_relations[-idx,];
idx<-which(inds_relations[,1]>inds_relations[,2]);
inds_relations[idx,]<-inds_relations[idx,c(2,1,3)];
idx<-which(inds_relations[,3]=="Second Order");
inds_relations[idx,3]<-"2";
idx<-which(inds_relations[,3]=="Third Order");
inds_relations[idx,3]<-"3";
idx<-which(inds_relations[,3]!="3"&inds_relations[,3]!="2");
inds_relations[idx,3]<-"1";
inds_relations<-unique(inds_relations);
write.table(inds_relations,file='inds_relations.filtered.2535.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE);
#######################################################################

##lack of first order and second order relatives, which are not included in the 2504 samples owing to population statistical requierment
#inds_relations<-read.table(file='inds_relations.filtered.2504.txt');
inds_relations<-read.table(file='inds_relations.filtered.2535.txt');
inds_relations[,1]<-as.character(inds_relations[,1]);
inds_relations[,2]<-as.character(inds_relations[,2]);
t2<-matrix(0,dim(inds_relations)[1],2);
t2[,1]<-match(inds_relations[,1],inds[,1]);
t2[,2]<-match(inds_relations[,2],inds[,1]);
inds_similarity[t2]<-inds_relations[,3];
inds_similarity[t2[,c(2,1)]]<-inds_relations[,3];
inds_similarity[cbind(1:ind_num,1:ind_num)]<-0;

write.table(inds_similarity,file='inds_similarity_matrix.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE);
write.table(inds_similarity,file='inds_similarity_matrix.incl.name.txt',sep='\t',eol='\n',append=FALSE,quote=FALSE,col.names=TRUE,row.names=TRUE);

inds_ranking_nums<-matrix(0,2504,10);
for(i in 1:2504){for(j in 1:10){inds_ranking_nums[i,j]<-length(which(inds_similarity[i,]==j-1));}}
inds_ranking_ceil<-matrix(0,2504,10);
inds_ranking_floor<-matrix(0,2504,10);
for(i in 1:2504){for(j in 1:10){inds_ranking_ceil[i,j]<-sum(inds_ranking_nums[i,1:j]);}}
for(i in 1:2504){for(j in 2:10){inds_ranking_floor[i,j]<-sum(inds_ranking_nums[i,1:(j-1)]);}}

save(inds_ranking_nums,inds_ranking_floor,inds_ranking_ceil,inds_similarity,file='inds_similarity_matrix.RData')

#calculate data rank
#in the final data_rank, zero means the distance between the individuals is in the proper ranking range, 
#negative value means how far the individual before the floor of proper ranking range (in terms of the relations between the two individuals)
#positive value means how far the individual after the ceil of proper ranking range (in terms of the relations between the two individuals)

delta_rank<-matrix(0,2504,2504);
for(i in 1:2504){delta_rank[i,]<-rank(delta[i,]);}
for(i in 1:2504){
	for(j in 1:2504){
		if(delta_rank[i,j]>inds_ranking_floor[i,inds_similarity[i,j]+1]&&delta_rank[i,j]<=inds_ranking_ceil[i,inds_similarity[i,j]+1]){
			delta_rank[i,j]<-0;
		}else if(delta_rank[i,j]<=inds_ranking_floor[i,inds_similarity[i,j]+1]){
			delta_rank[i,j]<-delta_rank[i,j]-inds_ranking_floor[i,inds_similarity[i,j]+1]-1;
		}else if(delta_rank[i,j]>inds_ranking_ceil[i,inds_similarity[i,j]+1]){
			delta_rank[i,j]<-delta_rank[i,j]-inds_ranking_ceil[i,inds_similarity[i,j]+1];
		}
	}
}

##a better implements
ind_num<-2504;
delta_rank<-t(apply(delta,1,rank));
t2<-cbind(rep(1:ind_num,rep(ind_num,ind_num)),rep(1:ind_num,ind_num));
delta_rank_floor<-matrix(0,ind_num,ind_num);
delta_rank_floor[t]<-delta_rank[t]-inds_ranking_floor[cbind(t[,1],inds_similarity[t]+1)]-1;
delta_rank_floor[delta_rank_floor>0]<-0;
delta_rank_ceil<-matrix(0,2504,2504);
delta_rank_ceil[t]<-delta_rank[t]-inds_ranking_ceil[cbind(t[,1],inds_similarity[t]+1)];
delta_rank_ceil[delta_rank_ceil<0]<-0;

delta_rank<-delta_rank_floor+delta_rank_ceil;

results<-matrix(0,10,4);
total_mean_offset<-sum(abs(delta_rank))/length(delta_rank);
for(i in 1:10){
	results[i,1]<-length(which(delta_rank==0&inds_similarity==i-1));
	results[i,2]<-length(which(inds_similarity==i-1));
	results[i,3]<-results[i,1]/results[i,2];
	results[i,4]<-sum(abs(delta_rank[inds_similarity==i-1]))/length(which(inds_similarity==i-1&delta_rank!=0));
}

delta_rank_random<-matrix(0,2504,2504);
for(i in 1:2504){delta_rank_random[i,]<-sample(1:2504,2504);}

#calculate random rank expectation

stats<-matrix(0,2535,10);
expectation<-matrix(0,2535,10);
for(i in 1:2535){
	stat<-table(inds_similarity[i,]);
	idx<-as.numeric(names(stat));
	stats[i,idx+1]<-stat;
	for(j in 1:10){
		expectation[i,j]<-(sum(stats[i,0:(j-1)])^2+(2535-sum(stats[i,0:j]))^2)/5070;
	}
}

expectation_col<-colSums(expectation*stats)/pmax(colSums(stats),1);
expectation_total<-mean(rowSums(expectation*stats)/2535);

#calculate random rank expectation for any number of individuals

#Avg. Random Misranked Positions AMRP
a<-read.table(file='../Desktop/sampling_plans_sample_num_map.txt',sep='\t',header=FALSE);
a[,1]<-as.character(a[,1]);
a[,2]<-as.character(a[,2]);
for(i in 3:12){
	a[,i]<-as.numeric(as.character(a[,i]));
}
popres<-read.table(file='../Desktop/pop_relations.txt',sep='\t',header=TRUE);
popres[,1]<-as.character(popres[,1]);
for(i in 2:27){
	popres[,i]<-as.numeric(as.character(popres[,i]));
}

expectation_total<-numeric(10);
expectation_col<-matrix(0,10,10);
for(k in 3:12){
	n<-a[27,k];
	stats<-matrix(0,n,10);
	expectation<-matrix(0,n,10);
	idx<-0;
	for(g in 1:26){
		for(i in 1:a[g,k]){
			idx<-idx+1;
			stats[idx,1]<-1;
			for(g2 in 1:26){
				if(popres[g,g2+1]==5){
					stats[idx,6]<-a[a[,2]==popres[g2,1],k]-1;
				}else{
					stats[idx,popres[g,g2+1]+1]<-stats[idx,popres[g,g2+1]+1]+a[a[,2]==popres[g2,1],k];
				}
			}
			for(j in 1:10){
				expectation[idx,j]<-(sum(stats[idx,0:(j-1)])^2+(n-sum(stats[idx,0:j]))^2)/(n*2);
			}
		}
	}
	expectation_total[k-2]<-mean(rowSums(expectation*stats)/n);
	expectation_col[k-2,]<-colSums(expectation*stats)/pmax(colSums(stats),1);
}
