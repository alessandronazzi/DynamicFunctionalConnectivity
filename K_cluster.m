
function [idx,centroids,upts] = K_cluster(upt_vmat,num_cluster)
     
    upt = cell2mat(upt_vmat);
    uptr = reshape(upt,[length(upt_vmat{1}),length(upt_vmat),size(upt_vmat,2)]);
    upts = transpose(reshape(uptr,[length(upt_vmat{1}),length(upt_vmat)*size(upt_vmat,2)]));
    
    [idx, centroids] = kmeans(upts,num_cluster);

end    
