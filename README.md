# LDA_genome

Source Codes:

generate_data.pl: Extract KGP data to variant information and individual genotypes.

merge_relate_data.pl: Merge 31 relatives data into the standard 2,504 dataset.

planning_sampling.R: Sampling variants for LDA training based on variant information.

sampling_and_insert_tab.pl: Extract genotypes data for the sampled variants.

calculate_inds_similarity_ref.R: Calculate reference ranking for KGP individuals.

modelling.R: LDA modelling.

calculate_pca_distance.R: PCA modelling.

calculate_hamming_distance.R: Hamming Distance modelling.

clustering_ranking_analysis.R: Ranking Score calculation.

plot.mds.Rplot.mds.R: MDS analysis and visualization.

Data access:

The 1000 genomes project dataset is downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

The 31 related individuals dataset is downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/related_samples_vcf/

High resolution figures: Figure 1-3, Figure S1-S3.
