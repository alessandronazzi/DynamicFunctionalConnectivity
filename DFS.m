

function [centroids] = DFS(timeseries_all, width_w, n_sbj)

    ex_check = exist('n_sbj','var');
    if ex_check == 0
        n_sbj = size(timeseries_all,3);
    else
    
    top_w = 1; 
    bot_w = width_w;
    n_window = size(timeseries_all,2) - bot_w + 1; 
    
    Mcorrg = []; %preallocate zeros to decrease runtime
    Mcorrs = [];
    
    % SLIDING-WINDOWS %
    
    for c=1:n_sbj
        SBJn = timeseries_all(:,:,c);
        for i=1:n_window
            timewindow_corr = corr(SBJn(:,top_w:bot_w)');
            Mcorrs = [Mcorrs; {timewindow_corr}];
            top_w = top_w+1;
            bot_w = bot_w+1;
        end    
        Mcorrg = [Mcorrg Mcorrs];
        Mcorrs = [];
        top_w = 1;
        bot_w = 86;
    end
    
    % EIGENDECOMPOSITION AND CORRESPONDING MATRICES %
    
    eigv_mat = [];
    mat_dim = length(Mcorrg) * n_sbj;
    
    for i=1:mat_dim
       [autvet, ~, ~] = svd(Mcorrg{i});
       eigv_mat = [eigv_mat; {autvet(:,1) * transpose(autvet(:,1))}];
    end
    
    eigv_mat = reshape(eigv_mat,length(Mcorrg),n_sbj);
    
    % UPPER-TRIANGULAR VECTORIZATION %
    
    upt_vmat = [];
    
    for c=1:mat_dim
        mat = eigv_mat{c};
        upt_v = mat(triu(true(size(mat)),1));
        upt_vmat = [upt_vmat; {upt_v}];
    end
    
    upt_vmat = reshape(upt_vmat,length(Mcorrg),n_sbj);
     
    % K-MEANS CLUSTERING %  
    
    upt = cell2mat(upt_vmat);
    uptr = reshape(upt,[length(upt_v),n_window,n_sbj]);
    upts = reshape(uptr,[length(upt_v),n_window*n_sbj]);
    
    [idx, centroids] = kmeans(upts',5);
    
    %  DFCs %
    
    DFCs = [];
    C = zeros(200);
    cent_cell = mat2cell(centroids, [1 1 1 1 1]); % cleener way to divide matrix in cells
    
    for i=1:length(cent_cell)
        c = centroids(i,:);
        C(triu(true(200),1)) = c;
        C = C+C';
        DFCs = [DFCs, {C}];
    end
    
    DFS1 = DFCs{1};
    DFS2 = DFCs{2};
    DFS3 = DFCs{3};
    DFS4 = DFCs{4};
    DFS5 = DFCs{5};

    figure()
    imagesc(DFS1)
    figure()
    imagesc(DFS2)
    figure()
    imagesc(DFS3)
    figure()
    imagesc(DFS4)
    figure()
    imagesc(DFS5)
    
    end
end









 