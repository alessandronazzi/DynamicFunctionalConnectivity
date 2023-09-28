%atlas = 'GordonLaumann';
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


if strcmp(atlas,'old')
    R=load([resultdir,'K',num2str(K),'.mat']);
    homotopic=R.all.metrics(:,1);
    dandmn=R.all.metrics(:,3);
    modularity=R.all.metrics(:,5);

else






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

if(strcmp(atlas,'GordonLaumann') || strcmp(atlas,'TianS1') || strcmp(atlas,'TianS4')) 
   load('/home/allegra/hcp_padova/results/GordonLaumann/NET_info_reduction2.mat');
end


B=NET_new.Index_reduction;

LR = zeros(90,1);

for i=1:37
   w=B{i};
   LR(w)=NET_new.LEFT_RIGHT(i); 
end

left = find(LR==1);
right = find(LR==2);

centroids = R.results.centroids;

centroid_mats = zeros(90,90,5);

idx = find(triu(ones(90),1));

homotopic = zeros(5,1);
dandmn = zeros(5,1);
modularity = zeros(5,1);

for i=1:5
    M=zeros(90);
    M(idx)=centroids(i,:);
    M=M+M';
    centroid_mats(:,:,i)=M;
    
    for m=1:10
        vl=find(NET_new.net_reduced==i & LR==1);
        vr=find(NET_new.net_reduced==i & LR==2);
    
        homotopic(i)=homotopic(i)+sum(sum(M(vl,vr)));
    end  
  
    danl = find(NET_new.net_reduced==6 & LR==1);
    danr = find(NET_new.net_reduced==6 & LR==1);
    dmnl = find(NET_new.net_reduced==8 & LR==2);
    dmnr = find(NET_new.net_reduced==8 & LR==2);
    
    dandmn(i)=sum(sum(M(danl,dmnl)));
    dandmn(i)=sum(sum(M(danr,dmnr)));
   
  
    
    W0 = M.*(M>0);

    m0 = Modul(NET_new.net_reduced,W0);    
    s0=sum(sum(W0));
    %B0 = W0 - 1./s0*sum(W0,2)*sum(W0,1);

    W1 = M.*(M<0);
    m1 = Modul(NET_new.net_reduced,-W1);
    s1=-sum(sum(W1));
    
    %B1 = W1 - 1./s1*sum(W1,2)*sum(W1,1);
    %B = B0/(s0+s1)-B1/(s0+s1);
    
    modularity(i) = (m0*s0-m1*s1)/(s0+s1);
    
close   
end
end


ax0=figure

ax1=subplot(1,3,1)
plot(homotopic)
%pbaspect([1 1 1]);
xlim([0.5,5.5])
title('homotopic FC')

ax2=subplot(1,3,2)
plot(dandmn)
%pbaspect([1 1 1]);
xlim([0.5,5.5])
title('dan-dmn FC')

ax3=subplot(1,3,3)
plot(modularity)
%pbaspect([1 1 1]);
xlim([0.5,5.5])
title('modularity')

pbaspect(ax1,[1 1 1])
pbaspect(ax2,[1 1 1])
pbaspect(ax3,[1 1 1])


saveas(gcf,['/home/allegra/hcp_padova/results/clean/graphs/cluster_metrics_',atlas,'.png'])

