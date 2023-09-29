function [optim_k,idx_opt,centroids_opt,mean_sil,std_sil] = optimalclust(upt_vmat,Npoints,Nrand,K)

   disp('Clustering');

    idx = {};

    for i=2:K
        disp(sprintf('k=%d',i));
        idx{i} = [];
        [idxs,centroids,upts] = K_cluster(upt_vmat,i);
        idx{i} = idxs;
        centroid{i} = centroids;
    end

    sil = {};

    disp('Compute silhouette');

    %tic
    for i=2:K
        disp(sprintf('k=%d',i));
        silh_vals=zeros(Nrand,1);
        indices = idx{i};
        parfor k=1:Nrand
            points=randperm(size(upts,1));
            points=points(1:Npoints);
            silh_vals(k)=mean(silhouette(upts(points,:),indices(points)));
        end
            sil{i}=silh_vals;
    end
    %toc

    for i=1:K
        mean_sil(i)=mean(sil{i},'omitnan');
        std_sil(i)=std(sil{i},'omitnan');
    end

    [~,optim_k] = max(mean_sil);
    idx_opt = idx{optim_k};
    centroids_opt = centroid{optim_k};

    %errorbar(mean_sil,std_sil);

end
