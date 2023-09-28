atlas1 = 'GordonLaumann';
%atlas1 = 'Schaefer';
%atlas1 = 'TianS1';
%atlas1 = 'TianS4';


atlas2 = 'GordonLaumann';
%atlas2 = 'Schaefer';
%atlas2 = 'TianS1';
%atlas2 = 'TianS4';
%atlas2 = 'old';

K1 = 5;  
do_pr1 = false;
glob1=true;
%distance1='Euclidean';
distance1='correlation';
do_downsample1=true;

K2 = 5;  
do_pr2 = true;
glob2=true;
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


nregions=90;
    
DFSs1 = [];
if(~strcmp(atlas1,'old'))
    centroids1=R1.results.centroids;

    for n=1:size(centroids1,1)
        C = zeros(nregions);
        c = centroids1(n,:);
        C(triu(true(nregions),1)) = c;
        C=C+C';
        DFSs1 = [DFSs1, {C}];
    end
else 
    for i=1:K1
        eval(['C=R1.all.DFS',num2str(i),';']);
        DFSs1 = [DFSs1, {C}];
    end
end

DFSs2 = [];
if(~strcmp(atlas2,'old'))
    centroids2=R2.results.centroids;

    for n=1:size(centroids2,1)
        C = zeros(nregions);
        c = centroids2(n,:);
        C(triu(true(nregions),1)) = c;
        C=C+C';
        DFSs2 = [DFSs2, {C}];
    end
else 
    for i=1:K1
        eval(['C=R2.all.DFS',num2str(i),';']);
        DFSs2 = [DFSs2, {C}];
    end
end


Conf=zeros(K1,K2);


for i=1:K1
    d1=DFSs1{i};
    for j=1:K2
	d2=DFSs2{j};
	idx=find(triu(ones(nregions),1));

        Conf(i,j)=corr(d1(idx),d2(idx));	

   end
end       


if do_downsample1 && strcmp(distance1,'correlation') && glob1 && ~do_pr1
    ylab='replication';
elseif do_downsample1 && strcmp(distance1,'correlation') && glob1 && do_pr1
    ylab='phase randomized';
elseif do_downsample1 && strcmp(distance1,'correlation') && ~glob1 && ~do_pr1
    ylab='no GSR';
elseif do_downsample1 && ~strcmp(distance1,'correlation') && glob1 && ~do_pr1
    ylab='Euclidean distance';
elseif ~do_downsample1 && strcmp(distance1,'correlation') && glob1 && ~do_pr1
    ylab='no downsampling';
end

if strcmp(atlas2,'old')
    xlab='original';
elseif do_downsample2 && strcmp(distance2,'correlation') && glob2 && ~do_pr2
    xlab='replication';
elseif do_downsample2 && strcmp(distance2,'correlation') && glob2 && do_pr2
    xlab='phase randomized';
elseif do_downsample2 && strcmp(distance2,'correlation') && ~glob2 && ~do_pr2
    xlab='no GSR';
elseif do_downsample2 && ~strcmp(distance2,'correlation') && glob2 && ~do_pr2
    xlab='Euclidean distance';
elseif ~do_downsample2 && strcmp(distance2,'correlation') && glob2 && ~do_pr2
    xlab='no downsampling';
end

heatmap(Conf)
xlabel(xlab);
ylabel(ylab);

saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/confusion_centroids_',ylab,'_',xlab,'.png'])


