function main

    function [upt_vmat] = projLEigv(timeseries_all,width_w,n_sbj)
    
        ex_check = exist('n_sbj','var');
        if ex_check == 0
            n_sbj = size(timeseries_all,3);
        else
        
            top_w = 1; 
            bot_w = width_w;
            n_window = size(timeseries_all,2) - bot_w + 1; 
            
            Mcorrg = []; 
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
        end 
    end    
    
    function [mean_sil,std_sil,optim_k] = optimalclust(upt_vmat,Npoints,Nrand,Nclustering)
    
        idx = {};
        
        for i=1:Nclustering
            idx{i} = [];
            [idxs,~,upts] = K_cluster(upt_vmat,i);
            idx{i} = idxs;
        end    
        
        sil = {};
        
        for i=2:Nclustering
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
    
        errorbar(mean_sil,std_sil);
    
        for i=1:length(mean_sil)
            if mean_sil(i) == max(mean_sil)
                optim_k = i;
            else
                continue
            end
        end    
    end    
    
    function [idx,centroids,upts] = K_cluster(upt_vmat,otpim_k)
         
        upt = cell2mat(upt_vmat);
        uptr = reshape(upt,[length(upt_vmat{1}),length(upt_vmat),size(upt_vmat,2)]);
        upts = transpose(reshape(uptr,[length(upt_vmat{1}),length(upt_vmat)*size(upt_vmat,2)]));
        
        [idx, centroids] = kmeans(upts,otpim_k);
    
    end    
    
    
    function [FC_dynamics1,FC_dynamics2,Fraction_times,Dwell_times,Transition_prob] = FC_dyn(idx)
    
        % FRACTION TIMES %
        
        [C,~,ic] = unique(idx);
        idx_counts = accumarray(ic,1);
        num_clusters = length(C);
    
        Fraction_times = zeros([5 1]);
        
        for i=1:length(idx_counts)
            Fraction_times(i) = idx_counts(i)/length(idx)*100;
        end 
    
        % DWELL TIMES %
        
        ideal_dwell = zeros([1 length(idx)]);
        DFS_sequence = [];
        seq_num = 1;
        
        for i=1:length(idx)
            if i == length(idx)
                ideal_dwell(seq_num) = ideal_dwell(seq_num)+1;
                if idx(end) == idx(end-1)
                    DFS_sequence = [DFS_sequence, ideal_dwell(seq_num)];
                else
                    seq_num = seq_num+1;
                end
            else
                ideal_dwell(seq_num) = ideal_dwell(seq_num)+1;
                if idx(i) == idx(i+1) 
                    continue
                else
                    DFS_sequence = [DFS_sequence, ideal_dwell(seq_num)];
                    seq_num = seq_num+1;
                end
            end
        end   
        
        idx_plus = [idx', idx(end)-1];
        diffidx = diff(idx_plus);
        seq_idx = idx(diffidx~=0);
        
        DFS_sequence = DFS_sequence';
        [seq_sorted, I] = sort(seq_idx);
        table_seq = [seq_sorted, DFS_sequence(I)];
    
        Dwell_times = [zeros(num_clusters,1)];
        seq_num = 1;
        c = 0;
    
        for i=1:length(table_seq)
             Dwell_times(seq_num) = Dwell_times(seq_num) + table_seq(i,2);
             c = c+1;
            if i == length(table_seq)           
                 Dwell_times(seq_num) = Dwell_times(seq_num)/c;
            else
                if table_seq(i) == table_seq(i+1)
                    continue            
                else
                    Dwell_times(seq_num) = Dwell_times(seq_num)/c;
                    c = 0;
                    seq_num = seq_num+1;
                end
            end    
        end
    %%
        % TRANSITION PROBABILITY %
        
        transitions = [zeros(1,num_clusters*num_clusters)];
        seq_num = 1;
        
        for i=1:num_clusters
            for j=1:num_clusters
                 if i == j
                     seq_num = seq_num + 1;
                     continue
                 else    
                     for k=1:length(seq_idx)
                         if k == length(seq_idx)
                             seq_num = seq_num + 1;
                             continue
                         else    
                            if seq_idx(k) == i && seq_idx(k+1) == j
                                transitions(seq_num) = transitions(seq_num) + 1;
                                continue
                            end
                         end   
                     end
                 end    
            end
        end    
        
        
        Transition_prob = [];
            
        for i=1:length(transitions)
            Transition_prob = [Transition_prob, (transitions(i)/(length(seq_idx)-1))*100];
        end
        
        Transition_prob = transpose(reshape(Transition_prob,[5 5]));
        
        % DISPLAY RESULTS %
        
        State = {'DFS1';'DFS2';'DFS3';'DFS4';'DFS5'};
        
        FC_dynamics1 = table(State,Fraction_times,Dwell_times);
        FC_dynamics2 = table(Transition_prob,'RowNames',State);
    end
end
        