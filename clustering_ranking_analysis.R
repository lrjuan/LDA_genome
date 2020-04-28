#args<-commandArgs(TRUE);
#files<-list.files(path=".",pattern=paste("sample",args[1],"*RData$",sep='.'));
files<-list.files(path=".",pattern="sample.100.9.*.RData$");
filenames<-strsplit(files,split='\\.');

library(slam);
library(topicmodels);

load('../inds_similarity_matrix.RData');
inds_similarity_all<-inds_similarity;

analysis_result<-list();
ranking_result<-list();

#kmeans_times<-10;
#k_range<-4:32;
kmeans_times<-5;
k_range<-26:26;

a_idx<-1;
r_idx<-1;

for (f in 1:length(files)) {
	load(files[f]);

#	inds<-inds_less;
#	ind_idx<-ind_idx_less;
	ind_num<-ind_num_less;
#	delta<-delta[ind_idx_less,ind_idx_less];
	inds_similarity<-inds_similarity_all[ind_idx_less,ind_idx_less];
	inds_ranking_nums<-matrix(0,ind_num_less,10);
	inds_ranking_ceil<-matrix(0,ind_num_less,10);
	inds_ranking_floor<-matrix(0,ind_num_less,10);
	for(i in 1:ind_num_less){for(j in 1:10){inds_ranking_nums[i,j]<-length(which(inds_similarity[i,]==j-1));}}
	for(i in 1:ind_num_less){for(j in 1:10){inds_ranking_ceil[i,j]<-sum(inds_ranking_nums[i,1:j]);}}
	for(i in 1:ind_num_less){for(j in 2:10){inds_ranking_floor[i,j]<-sum(inds_ranking_nums[i,1:(j-1)]);}}
	t2<-cbind(rep(1:ind_num_less,rep(ind_num_less,ind_num_less)),rep(1:ind_num_less,ind_num_less));

	delta_rank<-t(apply(delta,1,rank));
	delta_rank_floor<-matrix(0,ind_num,ind_num);
	delta_rank_floor[t2]<-delta_rank[t2]-inds_ranking_floor[cbind(t2[,1],inds_similarity[t2]+1)]-1;
	delta_rank_floor[delta_rank_floor>0]<-0;
	delta_rank_ceil<-matrix(0,ind_num,ind_num);
	delta_rank_ceil[t2]<-delta_rank[t2]-inds_ranking_ceil[cbind(t2[,1],inds_similarity[t2]+1)];
	delta_rank_ceil[delta_rank_ceil<0]<-0;
	delta_rank<-delta_rank_floor+delta_rank_ceil;

	delta_rank_random<-matrix(0,ind_num,ind_num);
	for(i in 1:ind_num){delta_rank_random[i,]<-sample(1:ind_num,ind_num);}
	delta_rank_floor<-matrix(0,ind_num,ind_num);
	delta_rank_floor[t2]<-delta_rank_random[t2]-inds_ranking_floor[cbind(t2[,1],inds_similarity[t2]+1)]-1;
	delta_rank_floor[delta_rank_floor>0]<-0;
	delta_rank_ceil<-matrix(0,ind_num,ind_num);
	delta_rank_ceil[t2]<-delta_rank_random[t2]-inds_ranking_ceil[cbind(t2[,1],inds_similarity[t2]+1)];
	delta_rank_ceil[delta_rank_ceil<0]<-0;
	delta_rank_random<-delta_rank_floor+delta_rank_ceil;

	for(i in c(1:3,5:9)){
		ranking_result[[r_idx]]<-c(as.numeric(filenames[[f]][2]),as.numeric(filenames[[f]][3]),i,length(which(delta_rank==0&inds_similarity==i))
							   ,length(which(delta_rank_random==0&inds_similarity==i))
							   ,length(which(inds_similarity==i))
							   ,sum(abs(delta_rank[inds_similarity==i]))/length(which(inds_similarity==i&delta_rank!=0))
							   ,sum(abs(delta_rank_random[inds_similarity==i]))/length(which(inds_similarity==i&delta_rank_random!=0)));
		r_idx<-r_idx+1;
	}

	inds_assignments<-list();
	inds_assignments[[1]]<-topics(lda);

	delta<-as.dist(delta);
	ks<-rep(k_range,rep(kmeans_times,length(k_range)));
	for(k_idx in 1:length(ks)){
		kc<-kmeans(delta,ks[k_idx]);
		inds_assignments[[k_idx+1]]<-kc$cluster;
	}
	ks<-c(0,ks);
	for(a in 1:length(inds_assignments)){
		inds_assignment<-numeric(2535);
		inds_assignment[ind_idx_less]<-inds_assignments[[a]];
		cluster_names<-unique(inds_assignments[[a]]);
		clusters<-list();
		for(i in 1:length(cluster_names)){clusters[[i]]<-which(inds_assignment==cluster_names[i]);}
		pops<-inds_pop;
	
		t<-t(combn(1:ind_num,2));
	
		P<-matrix(0,length(clusters),length(pops))
		R<-matrix(0,length(clusters),length(pops))
		F_measure<-matrix(0,length(clusters),length(pops))
		for(i in 1:length(clusters)){
			for(j in 1:length(pops)){
				P[i,j]<-length(intersect(pops[[j]],clusters[[i]]))/length(clusters[[i]]);
				R[i,j]<-length(intersect(pops[[j]],clusters[[i]]))/length(pops[[j]]);
				if(P[i,j]!=0 || R[i,j]!=0){
					F_measure[i,j]<-2*P[i,j]*R[i,j]/(P[i,j]+R[i,j]);
				}
			}
		}
		F_pops<-apply(F_measure,2,max);
		pops_weight<-sapply(pops,length)/ind_num;
	
		class_F<-sum(F_pops*pops_weight);  

		analysis_result[[a_idx]]<-c(as.numeric(filenames[[f]][2]),as.numeric(filenames[[f]][3]),ks[a],class_F);
		a_idx<-a_idx+1;
	}
}

topic_matrix<-t(as.data.frame(analysis_result[seq(1,length(analysis_result),length(ks))]));
rownames(topic_matrix)<-NULL;

analysis_result<-t(as.data.frame(analysis_result[-seq(1,length(analysis_result),length(ks))]));
rownames(analysis_result)<-NULL;
for(i in 2:kmeans_times){
	analysis_result[seq(1,dim(analysis_result)[1],kmeans_times),]<-analysis_result[seq(1,dim(analysis_result)[1],kmeans_times),]+analysis_result[seq(i,dim(analysis_result)[1],kmeans_times),];
}
analysis_result<-analysis_result[seq(1,dim(analysis_result)[1],kmeans_times),]/kmeans_times;
ranking_result<-t(as.data.frame(ranking_result));
rownames(ranking_result)<-NULL;

ranking_result_summary<-cbind(ranking_result[seq(1,dim(ranking_result)[1],8),1:2],matrix(0,dim(ranking_result)[1]/8,10));
for(i in 1:8){
	ranking_result_summary[,i+2]<-ranking_result[seq(i,dim(ranking_result)[1],8),7]*(ranking_result[seq(i,dim(ranking_result)[1],8),6]-ranking_result[seq(i,dim(ranking_result)[1],8),4])/ranking_result[seq(i,dim(ranking_result)[1],8),6];
	ranking_result_summary[,11]<-ranking_result_summary[,11]+ranking_result[seq(i,dim(ranking_result)[1],8),7]*(ranking_result[seq(i,dim(ranking_result)[1],8),6]-ranking_result[seq(i,dim(ranking_result)[1],8),4]);
	ranking_result_summary[,12]<-ranking_result_summary[,12]+ranking_result[seq(i,dim(ranking_result)[1],8),6];
}
ranking_result_summary[,11]<-ranking_result_summary[,11]/ranking_result_summary[,12];

ranking_random_summary<-cbind(ranking_result[seq(1,dim(ranking_result)[1],8),1:2],matrix(0,dim(ranking_result)[1]/8,10));
for(i in 1:8){
	ranking_random_summary[,i+2]<-ranking_result[seq(i,dim(ranking_result)[1],8),8]*(ranking_result[seq(i,dim(ranking_result)[1],8),6]-ranking_result[seq(i,dim(ranking_result)[1],8),5])/ranking_result[seq(i,dim(ranking_result)[1],8),6];
	ranking_random_summary[,11]<-ranking_random_summary[,11]+ranking_result[seq(i,dim(ranking_result)[1],8),8]*(ranking_result[seq(i,dim(ranking_result)[1],8),6]-ranking_result[seq(i,dim(ranking_result)[1],8),5]);
	ranking_random_summary[,12]<-ranking_random_summary[,12]+ranking_result[seq(i,dim(ranking_result)[1],8),6];
}
ranking_random_summary[,11]<-ranking_random_summary[,11]/ranking_random_summary[,12];

overall_results<-cbind(topic_matrix,analysis_result[,3:4],ranking_result_summary[,3:11],ranking_random_summary[,3:11]);

#save(topic_matrix,ranking_result,analysis_result,file=paste('analysis_result',args[1],'RData',sep='.'));
save(topic_matrix,ranking_result,analysis_result,overall_results,file=paste('analysis_result','RData',sep='.'));
	
q(save='no');
