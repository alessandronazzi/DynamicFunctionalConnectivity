
function [idx,idx_clean,centroids,upts] = K_cluster(upt_vmat,K,distance)
    upts = transpose(upt_vmat);
    
    if strcmp(distance,'Euclidean')
        [idx, centroids] = kmeans(upts,K,'MaxIter',1000);
    else
        [idx, centroids] = kmeans(upts,K,'Distance','correlation','MaxIter',1000);
    end

    idx_clean = idx;


    if strcmp(distance,'Euclidean')
        Dth=min(pdist(centroids));
    else
        Dth=min(pdist(centroids,'Correlation'));
    end



    for i=1:length(idx)
	X=[upts(i,:); centroids(idx(i),:)];
	D=pdist(X,'Correlation');
	if(D>Dth*1.00)
	    idx_clean(i)=0;
	end
    end

end
