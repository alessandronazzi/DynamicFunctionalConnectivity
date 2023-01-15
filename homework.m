%%

clc;
clearvars;

cd 'C:\Users\aless\OneDrive\Desktop\PNC\TIRO\Schaefer2018_200Parcels_17Networks_order'
load timeseries_all.mat

%% SLIDING-WINDOWS %%

n_sbj = 10; % = size(timeseries_all,3)
% width_w parametro funzione
top_w = 1; % = width_w - (width_w - 1)
bot_w = 86; % = width_w
n_window = size(timeseries_all,2) - bot_w; 

Mcorrg = []; %preallocate zeros to decrease runtime
Mcorrs = [];

for c=1:n_sbj
    SBJn = timeseries_all(:,:,c);
    for i=1:n_window
        timewindow_corr = corr(SBJn(:,top_w:bot_w)');
        Mcorrs = [Mcorrs; {timewindow_corr}];
        top_w = top_w+1;
        bot_w = bot_w+1;
    end    
    Mcorrg = [Mcorrg Mcorrs];
    Mcorrs = [];
    top_w = 1;
    bot_w = 86;
end

%% EIGENDECOMPOSITION AND CORRESPONDING MATRICES %%

eigv_mat = [];
mat_dim = length(Mcorrg) * n_sbj;

for i=1:mat_dim
    [autvet, eigval] = eig(Mcorrg{i});
    eigv_mat = [eigv_mat; {autvet(:,1) * transpose(autvet(:,1))}];
end

eigv_mat = reshape(eigv_mat,length(Mcorrg),n_sbj);

%% UPPER-TRIANGULAR VECTORIZATION %%

upt_vmat = [];

for c=1:mat_dim
    mat = eigv_mat{c};
    upt_v = mat(triu(true(size(mat))));
    upt_vmat = [upt_vmat; {upt_v}];
end

upt_vmat = reshape(upt_vmat,length(Mcorrg),n_sbj);
 
%% RESULTS OVERVIEW %%

for a = 1:20
    test_img = eigv_mat{a};
    imagesc(test_img);
    pause
end







 