# DynamicalFunctionalConnectivity
This repository contains the code employed for the Master's thesis titled "Dynamic connectivity patterns and cortico-subcortical interactions in the human brain" 

Language: MATLAB

Functions/: main functions and scripts

* Dfss.m: plots Dynamic Functional States for different atlases
* Fc_dyn.m: computes fraction times, dwell times and transition probabilities associated with the DFSs
* K_cluster.m: performs a K-means clustering with specified parameters
* Modul.m: computes classical Newman-Girvan modularity (credits E. le Martelot)
* behave.m: fits GLM with fraction and dwell times as predictors and different behavioral indices as dependent variables
* viol_dyn.m: plots fraction and dwell times distributions 
* cluster_metrics.m: computes and plots functional biomarkers associated with each DFSs
* compute_nmi.m: computes nmi scores for different conditions
* confm_cent.m: plots confusion matrices of the correlations between the centroids of different conditions
* confm_ind.m: plots confusion matrices of the cluster's assignations of different conditions
* main.m: retrieves the data and performs all the analysis steps calling different custom functions
* meandfc_eigv.m: plots one example of the average connectivity over time for cortical and subcortical networks (GordonLaumann cortical parcellation and FreeSurfer subcortical parcellation)
* nmi.m: computes normalized mutual information of two discrete variables (credits Mo Chen)
* optimalclust.m: computes silhouette values for different K, evaluating the optimal number of clusters for the K-means clustering
* phase_rand.m: performs phase randomization of the data
* plot_dyn.m: plots the distributions of fraction and dwell times across subjects
* preprocess.m: performs all the additional preprocessing steps included in the analysis
* probacdiff_eigv.m: plots the probability distribution of the absolute values of connectivity differences between consecutive sliding windows and the cumulative density function of the conditioned probability of subcortical connectivity reorganization, given a cortical connectivity reorganization.
* projLEigv.m: performs sliding-window temporal correlation returning the upper-triangular portion of the matrix reconstructed from the principal eigenvector of each sliding window
