
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
        















    

