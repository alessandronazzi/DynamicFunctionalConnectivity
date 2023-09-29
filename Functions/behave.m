atlas = 'GordonLaumann';
%atlas = 'Schaefer';
%atlas = 'TianS1';
%atlas = 'TianS4';

do_pr = false;
K=5;

resultdir=['/home/allegra/hcp_padova/results/clean/',atlas,'/'];

if do_pr
    R=load([resultdir,'K',num2str(K),'_pr.mat']);
else
    R=load([resultdir,'K',num2str(K),'.mat']);
end


%function [] = behave(idx)

bt = readtable('/home/allegra/hcp_padova/scripts/DynamicFunctionalConnectivity/beha.csv');



n_sbj = height(bt);
n_window = 746;
ss = zeros(1, n_sbj) + n_window;
idx_sbj = mat2cell(R.results.idx,ss,1);

ft_all = zeros(n_sbj,K);
dt_all = zeros(n_sbj,K);


for i=1:length(idx_sbj)
    sbj = idx_sbj{i};
    [ft,dt,~] = FC_dyn(sbj,K);
    ft_all(i,:)=ft;
    dt_all(i,:)=dt;
end


i1=find(strcmpi(bt.Properties.VariableNames,'SCPT_SEN')); %'MMSE_Score'));
i2=find(strcmpi(bt.Properties.VariableNames,'PSQI_Score'));
i3=find(strcmpi(bt.Properties.VariableNames,'PMAT24_A_CR'));
i4=find(strcmpi(bt.Properties.VariableNames,'ProcSpeed_Unadj'));
i5=find(strcmpi(bt.Properties.VariableNames,'MMSE_Score'));

used_vars = [i1 i2 i3 i4 i5];
nvars= length(used_vars);

btred = bt(:,used_vars);

Rsquare = zeros(nvars,1);
Pvals = zeros(nvars,1);

Nperm = 10000;

for n=1:nvars
    y=table2array(btred(:,n));
    used_sbj = find(~isnan(y(:,1)));
    y=y(used_sbj);
   
    X=ft_all(used_sbj,1:(K-1));
    
    [b,dev,stats]=glmfit(X,y);

    yp = X*b(2:end)+b(1);
    
    R= corrcoef(y,yp);
    Rsquare(n)=R(1,2)^2;
    %Pvals(n)=Pr(1,2);
    
    Rperm = zeros(Nperm,1);
    Rperm(1)=R(1,2);
    
    for r=2:Nperm
        
        ys = y(randperm(length(y))); 
        [bs,devs,statss]=glmfit(X,ys);

        yps = X*bs(2:end)+bs(1);
        Rs = corrcoef(ys,yps);
        Rperm(r)=Rs(1,2); 
    end
    
    Pvals(n) = length(find(Rperm>=Rperm(1)))/Nperm;
end
