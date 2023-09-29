atlas = 'GordonLaumann';
%%atlas = 'TianS1';
%atlas = 'TianS4';
%atlas= 'Schaefer';

K = 5;  
do_pr = false;
glob=false;
%distance='Euclidean';
distance='correlation';
do_downsample=true;


if strcmp(atlas,'GordonLaumann')
   datadir = '/home/allegra/hcp_padova/data/';
   files = dir([datadir,'*GLS.ptseries.mat']);
elseif strcmp(atlas,'Schaefer')
   datadir = '/home/allegra/hcp_padova/data/schaefer/';
   files = dir([datadir,'*.mat']);
elseif strcmp(atlas,'TianS1')
   datadir = '/home/allegra/hcp_padova/data/Tian_S1/';
   files = dir([datadir,'*.mat']);
else
   datadir = '/home/allegra/hcp_padova/data/Tian_S4/';
   files = dir([datadir,'*.mat']);   
end

%files = dir([datadir,'*.mat']);
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
else
    fname = [datadir,char(sbjs(1)),'.REST1.LR.GLTianS4.ptseries.mat'];
end
    
data_LR = load(fname);

[nregions,ntimepoints]=size(data_LR.tseries);
ntimepoints=2*ntimepoints;
nmax = 80;

timeseries_all = zeros(nregions,ntimepoints,nmax);

for n=1:nmax %nsbjs
    if strcmp(atlas,'GordonLaumann')
    	fname = [datadir,char(sbjs(n)),'.REST1.LR.GLS.ptseries.mat'];
    elseif strcmp(atlas,'Schaefer')
    	fname = [datadir,char(sbjs(n)),'.REST1.LR.SchaeferS.ptseries.mat'];
    elseif strcmp(atlas,'TianS1')
        fname = [datadir,char(sbjs(n)),'.REST1.LR.GLTianS1.ptseries.mat'];
    else
        fname = [datadir,char(sbjs(n)),'.REST1.LR.GLTianS4.ptseries.mat'];   
    end
    data_LR = load(fname);
    if strcmp(atlas,'GordonLaumann')
    	fname = [datadir,char(sbjs(n)),'.REST1.RL.GLS.ptseries.mat'];
    elseif strcmp(atlas,'Schaefer')
    	fname = [datadir,char(sbjs(n)),'.REST1.RL.SchaeferS.ptseries.mat'];
    elseif strcmp(atlas,'TianS1')
        fname = [datadir,char(sbjs(n)),'.REST1.RL.GLTianS1.ptseries.mat'];
    else
        fname = [datadir,char(sbjs(n)),'.REST1.RL.GLTianS4.ptseries.mat'];    
    end
    data_RL = load(fname);
    timeseries_all(:,:,n) = [data_LR.tseries(:,:), data_RL.tseries(:,:)];
end



%% PHASE RANDOMIZATION %%


if do_pr
    [timeseries_all] = phase_rand(timeseries_all);
end



%%

load('/home/allegra/hcp_padova/results/GordonLaumann/NET_info_reduction2.mat')
    

if glob
    [timeseries_all] = preprocess(timeseries_all,atlas);
else
    [timeseries_all] = preprocess(timeseries_all,atlas,'noglob');
end


% DOWNSAMPLING %

TR=0.71;
width_w = 84;

if(do_downsample)
    width_w = 28;  %floor(60/TR);

    w_timeseries = size(timeseries_all,1);
    l_timeseries = size(timeseries_all,2)/3;
    n_sbj = size(timeseries_all,3);

    timeseries_red = zeros(w_timeseries,l_timeseries,n_sbj);

    for i=1:size(timeseries_all,3)
        timeseries_red(:,:,i) = timeseries_all(:,1:3:l_timeseries*3,i);

    end    

    timeseries_all = timeseries_red;

end

%%% Compute principal eigenvector

[upt_vmat] = projLEigv(timeseries_all,width_w,nmax);
%[mat_vect] = projLEigv(timeseries_all,width_w,n_sbj);
sss 
Npoints = 100;
Nrand = 1000;


%[optim_k,idx_opt,centroids_opt,mean_silh,~] = optimalclust(upt_vmat,Npoints,Nrand,K);
[idx,centroids,upts] = K_cluster(upt_vmat,K,distance);


results.timeseries_all = timeseries_all;
results.centroids = centroids;
results.idx = idx;
%results.upts = upts(1:773*10,:);

[Fraction_times,Dwell_times,Transition_prob] = FC_dyn(idx,K);
    
results.fraction_times=Fraction_times;
results.dwell_times=Dwell_times;
results.Transition_prob=Transition_prob;

resultdir=['/home/allegra/hcp_padova/results/clean/',atlas,'/'];

if do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    save([resultdir,'K',num2str(K),'.mat'],'results');
elseif do_downsample && strcmp(distance,'correlation') && glob && do_pr
    save([resultdir,'K',num2str(K),'_pr.mat'],'results');
elseif do_downsample && strcmp(distance,'correlation') && ~glob && ~do_pr
    save([resultdir,'K',num2str(K),'_noglob.mat'],'results');
elseif do_downsample && ~strcmp(distance,'correlation') && glob && ~do_pr
    save([resultdir,'K',num2str(K),'_',distance,'.mat'],'results');
elseif ~do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    save([resultdir,'K',num2str(K),'_allpts.mat'],'results');
end


