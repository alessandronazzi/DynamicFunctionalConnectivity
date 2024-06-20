
function [optim_k,idx_opt,centroids_opt,mean_sil,std_sil] = optimalclust(upts,idx,centroid,Npoints,Nrand,K)

    sil = {};
    disp('Compute silhouette');
    
    for i=2:K
        disp(sprintf('k=%d',i));
        silh_vals=zeros(Nrand,1);
        indices = idx{i};
        parfor k=1:Nrand
            points=randperm(size(upts,1));
            points=points(1:Npoints);
            silh_vals(k)=mean(silhouette(upts(:,points)',indices(points)));
        end
	    sil{i}=silh_vals;
    end
    
    for i=1:K
        mean_sil(i)=mean(sil{i},'omitnan');
        std_sil(i)=std(sil{i},'omitnan');
    end
    
    [~,optim_k] = max(mean_sil);
    idx_opt = idx{optim_k};
    centroids_opt = centroid{optim_k};

end                             
