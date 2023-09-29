
function [idx,centroids,upts] = K_cluster(upt_vmat,K,distance)
    upts = transpose(upt_vmat);
    
    if strcmp(distance,'Euclidean')
        [idx, centroids] = kmeans(upts,K,'MaxIter',1000);
    else
        [idx, centroids] = kmeans(upts,K,'Distance','correlation','MaxIter',1000);
    end

end
