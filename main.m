datadir = '/home/allegra/hcp_padova/data/';
files = dir([datadir,'*.mat']);
nfiles = length(files);
sbjs = {};

for n=1:nfiles
   x=files(n).name;
   sbjs{n}=x(1:6);   
end

sbjs=unique(sbjs);
nsbjs=length(sbjs);

ntimepoints = 10*2;
nregions = 343;

nmax = 10;

timeseries_all = zeros(ntimepoints,nregions,nmax);

for n=1:nmax %nsbjs
    fname = [datadir,char(sbjs(n)),'.REST1.LR.GLS.ptseries.mat'];
    data_LR = load(fname);
    fname = [datadir,char(sbjs(n)),'.REST1.RL.GLS.ptseries.mat'];
    data_RL = load(fname);
    timeseries_all(:,:,n) = [data_LR.tseries(:,1:10),data_RL.tseries(:,1:10)]';
end


TR=0.71;
width_w =5; %floor(60/TR);

%%% Compute principal eigenvactors for sliding windows
[upt_vmat] = projLEigv(timeseries_all,width_w);    

Npoints = 100;
Nrand = 100;
K = 3;

[mean_sil,std_sil] = optimalclust(upt_vmat,Npoints,Nrand,K);

optim_k=2;

[idx,centroids,upts] = K_cluster(upt_vmat,optim_k);
 
[FC_dynamics1,FC_dynamics2,Fraction_times,Dwell_times,Transition_prob] = FC_dyn(idx);
    
disp(FC_dynamics1);
disp(FC_dynamics2);

