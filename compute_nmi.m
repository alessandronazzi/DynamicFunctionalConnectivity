%atlas1 = 'Schaefer';
%atlas1 = 'TianS1';
%atlas1 = 'TianS4';

atlas2 = 'GordonLaumann';
%atlas2 = 'Schaefer';
%atlas2 = 'TianS1';
%atlas2 = 'TianS4';

K1 = 5;  
do_pr1 = false;
glob1=true;
%distance1='Euclidean';
distance1='correlation';
do_downsample1=true;

K2 = 5;  
do_pr2 = false;
glob2=false; %true;
%distance2='Euclidean';
distance2='correlation';
do_downsample2=true;


resultdir1=['/home/allegra/hcp_padova/results/clean/',atlas1,'/'];
resultdir2=['/home/allegra/hcp_padova/results/clean/',atlas2,'/'];


if do_downsample1 && strcmp(distance1,'correlation') && glob1 && ~do_pr1
    R1=load([resultdir1,'K',num2str(K1),'.mat']);
elseif do_downsample1 && strcmp(distance1,'correlation') && glob1 && do_pr1
    R1=load([resultdir1,'K',num2str(K1),'_pr.mat']);
elseif do_downsample1 && strcmp(distance1,'correlation') && ~glob1 && ~do_pr1
    R1=load([resultdir1,'K',num2str(K1),'_noglob.mat']);
elseif do_downsample1 && ~strcmp(distance1,'correlation') && glob1 && ~do_pr1
    R1=load([resultdir1,'K',num2str(K1),'_',distance1,'.mat']);
elseif ~do_downsample1 && strcmp(distance1,'correlation') && glob1 && ~do_pr1
    R1=load([resultdir1,'K',num2str(K1),'_allpts.mat']);
end


if do_downsample2 && strcmp(distance2,'correlation') && glob2 && ~do_pr2
    R2=load([resultdir2,'K',num2str(K2),'.mat']);
elseif do_downsample2 && strcmp(distance2,'correlation') && glob2 && do_pr2
    R2=load([resultdir2,'K',num2str(K2),'_pr.mat']);
elseif do_downsample2 && strcmp(distance2,'correlation') && ~glob2 && ~do_pr2
    R2=load([resultdir2,'K',num2str(K2),'_noglob.mat']);
elseif do_downsample2 && ~strcmp(distance2,'correlation') && glob2 && ~do_pr2
    R2=load([resultdir2,'K',num2str(K2),'_',distance2,'.mat']);
elseif ~do_downsample2 && strcmp(distance2,'correlation') && glob2 && ~do_pr2
    R2=load([resultdir2,'K',num2str(K2),'_allpts.mat']);
end

nsbj1=size(R1.results.timeseries_all,3);
idxsub1=[];
for i=1:nsbj1
   idxsub1=[idxsub1 (2234*(i-1)+[1:3:1117 1118:3:2234]) ];
end

nsbj2=size(R2.results.timeseries_all,3);
idxsub2=[];
for i=1:nsbj2
   idxsub2=[idxsub2 (2234*(i-1)+[1:3:1117 1118:3:2234]) ];
end



if do_downsample1 && ~do_downsample2
	MI = nmi(R1.results.idx,R2.results.idx(idxsub2))
elseif ~do_downsample1 && do_downsample2
	MI = nmi(R1.results.idx(idxsub1),R2.results.idx)
elseif ~do_downsample1 && ~do_downsample2
	MI = nmi(R1.results.idx(idxsub1),R2.results.idx(idxsub2))
else do_downsample1 && do_downsample2
    MI = nmi(R1.results.idx,R2.results.idx)
end
