function [mean_sil,std_sil] = optimalclust(upt_vmat,Npoints,Nrand,Nclustering)

   disp('Clustering');

    idx = {};
    
    for i=1:Nclustering
        disp(sprintf('k=%d',i));
        idx{i} = [];
        [idxs,~,upts] = K_cluster(upt_vmat,i);
        idx{i} = idxs;
    end    
    
    sil = {};
    
    disp('Compute silhouette');
    for i=2:Nclustering
        disp(sprintf('k=%d',i));
        sil{i}=[];
        indices = idx{i};
        for k=1:Nrand
            points=randperm(size(upts,1));
            points=points(1:Npoints);
            sil{i} =[sil{i};mean(silhouette(upts(points,:),indices(points)))];
        end
	    
    end
    
    for i=1:Nclustering
        mean_sil(i)=mean(sil{i},'omitnan');
        std_sil(i)=std(sil{i},'omitnan');
    end

    %errorbar(mean_sil,std_sil);

end    



                                
