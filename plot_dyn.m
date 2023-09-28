%%%%%atlas = 'GordonLaumann';
%atlas = 'Schaefer';
%atlas = 'TianS1';
%atlas = 'TianS4';
atlas= 'old';

K = 5;
do_pr = false;
glob=true;
%distance='Euclidean';
distance='correlation';
do_downsample=true;


resultdir=['/home/allegra/hcp_padova/results/clean/',atlas,'/'];




if do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    R=load([resultdir,'K',num2str(K),'.mat']);
elseif do_downsample && strcmp(distance,'correlation') && glob && do_pr
    R=load([resultdir,'K',num2str(K),'_pr.mat']);
elseif do_downsample && strcmp(distance,'correlation') && ~glob && ~do_pr
    R=load([resultdir,'K',num2str(K),'_noglob.mat']);
elseif do_downsample && ~strcmp(distance,'correlation') && glob && ~do_pr
    R=load([resultdir,'K',num2str(K),'_',distance1,'.mat']);
elseif ~do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    R=load([resultdir,'K',num2str(K),'_allpts.mat']);
end


if do_downsample
  l=800;
  w=28;
else
  l=2400;
  w=84;
end

if ~strcmp(atlas,'old')
    [ft,dt]=boxp_dyn(80,R.results.idx,l,w);
else
    ft=R.all.fractiontimes;
    dt=R.all.fractiontimes;
end


figure()
boxplot(ft);
title('Fraction times')

if do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/ft_',atlas,'_K',num2str(K),'.png'])
elseif do_downsample && strcmp(distance,'correlation') && glob && do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/ft_',atlas,'_K',num2str(K),'_pr.png'])
elseif do_downsample && strcmp(distance,'correlation') && ~glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/ft_',atlas,'_K',num2str(K),'_noglob.png'])
elseif do_downsample && ~strcmp(distance,'correlation') && glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/ft_',atlas,'_K',num2str(K),'_Euclidean.png'])
elseif ~do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/ft_',atlas,'_K',num2str(K),'_allpts.png'])
end



figure()
boxplot(dt);
title('Dwell times')

if do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/dt_',atlas,'_K',num2str(K),'.png'])
elseif do_downsample && strcmp(distance,'correlation') && glob && do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/dt_',atlas,'_K',num2str(K),'_pr.png'])
elseif do_downsample && strcmp(distance,'correlation') && ~glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/dt_',atlas,'_K',num2str(K),'_noglob.png'])
elseif do_downsample && ~strcmp(distance,'correlation') && glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/dt_',atlas,'_K',num2str(K),'_Euclidean.png'])
elseif ~do_downsample && strcmp(distance,'correlation') && glob && ~do_pr
    saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/dt_',atlas,'_K',num2str(K),'_allpts.png'])
end
