function [DFSs] = Dfss(centroids)

    DFSs = [];
    C = eye(200);

    for i=1:size(centroids,1)
        c = centroids(i,:);
        C(triu(true(200),1)) = c;
        U = triu(C,1);
        C=C+U';
        DFSs = [DFSs, {C}];
    end
    
    for i=1:length(DFSs)
        name = "DFS"+i;
        figure()
        imagesc(DFSs{i})
        title(name)
    end    
end