atlas = 'GordonLaumann';
%atlas = 'TianS1';
%atlas = 'TianS4';
%atlas= 'Schaefer';
%atlas = 'Buckner'
%atlas= 'Morel';

K = 5;  
do_pr = false;
glob=true;
%distance='Euclidean';
distance='correlation';
do_downsample=false;

if strcmp(atlas,'GordonLaumann')
   datadir = '/home/allegra/hcp_padova/data/GLS/';
   files = dir([datadir,'*GLS.ptseries.mat']);
elseif strcmp(atlas,'Schaefer')
   datadir = '/home/allegra/hcp_padova/data/Schaefer/';
   files = dir([datadir,'*REST1*Schaefer*.mat']);
elseif strcmp(atlas,'TianS1')
   datadir = '/home/allegra/hcp_padova/data/Tian_S1/';
   files = dir([datadir,'*.mat']);
elseif strcmp(atlas,'TianS4')
   datadir = '/home/allegra/hcp_padova/data/Tian_S4/';
   files = dir([datadir,'*.mat']);   
elseif strcmp(atlas,'Buckner')
   datadir = '/home/allegra/hcp_padova/data/';
   files = dir([datadir,'*GLTianS1Buck.ptseries.mat']);
elseif strcmp(atlas,'Morel')
   datadir = '/home/allegra/hcp_padova/data/';
   files = dir([datadir,'*GLTianS1Morel.ptseries.mat']);
end

nfiles = length(files);
sbjs = {};

for n=1:nfiles
   x=files(n).name;
   sbjs{n}=x(1:6);   
end

sbjs=unique(sbjs);
n_sbj=length(sbjs);

if strcmp(atlas,'GordonLaumann')
    fname = [datadir,char(sbjs(1)),'.REST1.LR.GLS.ptseries.mat'];
elseif strcmp(atlas,'Schaefer')
    fname = [datadir,char(sbjs(1)),'.REST1.LR.SchaeferS.ptseries.mat'];
elseif strcmp(atlas,'TianS1')
    fname = [datadir,char(sbjs(1)),'.REST1.LR.GLTianS1.ptseries.mat'];
elseif strcmp(atlas,'TianS1')
    fname = [datadir,char(sbjs(1)),'.REST1.LR.GLTianS4.ptseries.mat'];
elseif strcmp(atlas,'Buckner')
    fname = [datadir,char(sbjs(1)),'.REST1.LR.GLTianS1Buck.ptseries.mat'];
elseif strcmp(atlas,'Morel')
    fname = [datadir,char(sbjs(1)),'.REST1.LR.GLTianS1Morel.ptseries.mat'];
end
    
data_LR = load(fname);
[nregions,ntimepoints]=size(data_LR.tseries);
ntimepoints=2*ntimepoints;
nmax = n_sbj;
timeseries_all = zeros(nregions,ntimepoints,nmax);
exc_sbj=[];

for n=1:nmax 
    if strcmp(atlas,'GordonLaumann')
    	fname = [datadir,char(sbjs(n)),'.REST1.LR.GLS.ptseries.mat'];
    elseif strcmp(atlas,'Schaefer')
    	fname = [datadir,char(sbjs(n)),'.REST1.LR.SchaeferS.ptseries.mat'];
    elseif strcmp(atlas,'TianS1')
        fname = [datadir,char(sbjs(n)),'.REST1.LR.GLTianS1.ptseries.mat'];
    elseif strcmp(atlas,'TianS4')
        fname = [datadir,char(sbjs(n)),'.REST1.LR.GLTianS4.ptseries.mat'];   
    elseif strcmp(atlas,'Buckner')
    	fname = [datadir,char(sbjs(n)),'.REST1.LR.GLTianS1Buck.ptseries.mat'];
    elseif strcmp(atlas,'Morel')
    	fname = [datadir,char(sbjs(n)),'.REST1.LR.GLTianS1Morel.ptseries.mat'];
    end
    data_LR = load(fname);
    if strcmp(atlas,'GordonLaumann')
    	fname = [datadir,char(sbjs(n)),'.REST1.RL.GLS.ptseries.mat'];
    elseif strcmp(atlas,'Schaefer')
    	fname = [datadir,char(sbjs(n)),'.REST1.RL.SchaeferS.ptseries.mat'];
    elseif strcmp(atlas,'TianS1')
        fname = [datadir,char(sbjs(n)),'.REST1.RL.GLTianS1.ptseries.mat'];
    elseif strcmp(atlas,'TianS4')
        fname = [datadir,char(sbjs(n)),'.REST1.RL.GLTianS4.ptseries.mat'];    
    elseif strcmp(atlas,'Buckner')
    	fname = [datadir,char(sbjs(n)),'.REST1.RL.GLTianS1Buck.ptseries.mat'];
    elseif strcmp(atlas,'Morel')
    	fname = [datadir,char(sbjs(n)),'.REST1.RL.GLTianS1Morel.ptseries.mat'];
    end
    
    data_RL = load(fname);
    
    if size(data_LR.tseries,2)<1200 || size(data_RL.tseries,2)<1200
        exc_sbj=[exc_sbj,n];
        continue
    else
        timeseries_all(:,:,n) = [data_LR.tseries(:,:), data_RL.tseries(:,:)];

    end

end

timeseries_all(:,:,exc_sbj)=[];
n_sbj=size(timeseries_all,3);
nmax=n_sbj;

%% PHASE RANDOMIZATION %%


if do_pr
    [timeseries_all] = phase_rand(timeseries_all);
end


%% GLOBAL SIGNAL REGRESSION %%


if glob
    [timeseries_all] = preprocess(timeseries_all,atlas);
else
    [timeseries_all] = preprocess(timeseries_all,atlas,'noglob');
end


%% DOWNSAMPLING %%

TR=0.71;
width_w = 84;
short='';
if width_w<84
	    short='_short';
end
if width_w<42
	    short='_short21';
end



if(do_downsample)
    width_w = width_w/3; 
    w_timeseries = size(timeseries_all,1);
    l_timeseries = size(timeseries_all,2)/3;
    n_sbj = size(timeseries_all,3);
    timeseries_red = zeros(w_timeseries,l_timeseries,n_sbj);
    for i=1:size(timeseries_all,3)
        timeseries_red(:,:,i) = timeseries_all(:,1:3:l_timeseries*3,i);

    end    
    timeseries_all = timeseries_red;
end

%% EIGENDECOMPOSITION %%

[upt_vmat,eigv_mat] = projLEigv(timeseries_all,width_w,n_sbj);

%% K-MEANS CLUSTERING %%

[idx,idx_clean,centroids,upts] = K_cluster(upt_vmat,K,distance);

%% SAVE RESULTS %%

results.centroids = centroids;
results.idx = idx;
resultdir=['/home/allegra/hcp_padova/results/final/',atlas,'/'];

if ~isdir(resultdir)
   mkdir(resultdir);
end

if do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'.mat'],'results');
elseif do_downsample && strcmp(distance,'correlation') && glob && do_pr 
      save([resultdir,'K',num2str(K),short,'_pr.mat'],'results');
elseif do_downsample && strcmp(distance,'correlation') && ~glob && ~do_pr  
      save([resultdir,'K',num2str(K),short,'_noglob.mat'],'results');
elseif do_downsample && ~strcmp(distance,'correlation') && glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'_',distance,'.mat'],'results');
elseif ~do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'_allpts.mat'],'results');
    w_timeseries = size(timeseries_all,1);
    l_timeseries = size(timeseries_all,2)/3;
    n_sbj = size(timeseries_all,3);
    timeseries_red = zeros(w_timeseries,l_timeseries,n_sbj);
    for i=1:size(timeseries_all,3)
        timeseries_red(:,:,i) = timeseries_all(:,1:3:l_timeseries*3,i);

    end
    timeseries_all = timeseries_red;
end

%% EIGENDECOMPOSITION %%

[upt_vmat,eigv_mat] = projLEigv(timeseries_all,width_w,n_sbj);

%% K-MEANS CLUSTERING %%

[idx,idx_clean,centroids,upts] = K_cluster(upt_vmat,K,distance);

%% SAVE RESULTS %%

results.centroids = centroids;
results.idx = idx;
resultdir=['/home/allegra/hcp_padova/results/final/',atlas,'/'];

if ~isdir(resultdir)
   mkdir(resultdir);
end

if do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'.mat'],'results');
elseif do_downsample && strcmp(distance,'correlation') && glob && do_pr
      save([resultdir,'K',num2str(K),short,'_pr.mat'],'results');
elseif do_downsample && strcmp(distance,'correlation') && ~glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'_noglob.mat'],'results');
elseif do_downsample && ~strcmp(distance,'correlation') && glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'_',distance,'.mat'],'results');
elseif ~do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
      save([resultdir,'K',num2str(K),short,'_allpts.mat'],'results');
end



