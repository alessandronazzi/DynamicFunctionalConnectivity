function [FC_dynamics1,FC_dynamics2] = FCdyn(idx)

    % FRACTION TIMES %
    
    [C,~,ic] = unique(idx);
    idx_counts = accumarray(ic,1);
    idx_occurence= [C, idx_counts];
    
    occurrence_prcntg = zeros([1 5]);
    
    for i=1:length(idx_counts)
        occurrence_prcntg(i) = idx_counts(i)/length(idx)*100;
    end 
    
    Fraction_times = occurrence_prcntg';
    
    % DWELL TIMES %
    
    ideal_dwell = zeros([1 length(idx)]);
    DFS_sequence = [];
    seq_num = 1;
    
    for i=1:length(idx)
        if i == length(idx)
            if idx(end) == idx(end-1)
                ideal_dwell(seq_num) = ideal_dwell(seq_num)+1;
                DFS_sequence = [DFS_sequence, ideal_dwell(seq_num)];
            else
                seq_num = seq_num+1;
                ideal_dwell(seq_num) = ideal_dwell(seq_num)+1;
            end
        else
            if idx(i) == idx(i+1) 
                ideal_dwell(seq_num) = ideal_dwell(seq_num)+1;
            else
                ideal_dwell(seq_num) = ideal_dwell(seq_num)+1;
                DFS_sequence = [DFS_sequence, ideal_dwell(seq_num)];
                seq_num = seq_num+1;
            end
        end
    end   
    
    idx_plus = [idx', idx(end)-1];
    diffidx = diff(idx_plus);
    seq_idx = idx(diffidx~=0);
    
    table_seq = [seq_idx, DFS_sequence'];
    
    occ_1 = [];
    occ_2 = [];
    occ_3 = [];
    occ_4 = [];
    occ_5 = [];
    
    
    for i=1:length(table_seq)
        if table_seq(i) == 1
           occ_1 = [occ_1, table_seq(i,2)];
        elseif table_seq(i) == 2
             occ_2 = [occ_2, table_seq(i,2)]; 
        elseif table_seq(i) == 3
             occ_3 = [occ_3, table_seq(i,2)];
        elseif table_seq(i) == 4
             occ_4 = [occ_4, table_seq(i,2)];
        elseif table_seq(i) == 5
             occ_5 = [occ_5, table_seq(i,2)];
        end
    end    
    
    Dwell_times = [mean(occ_1); mean(occ_2); mean(occ_3); mean(occ_4); mean(occ_5)];
    
    % TRANSITION PROBABILITY %
    
    trans_12 = 0;
    trans_21 = 0;
    trans_23 = 0;
    trans_32 = 0;
    trans_34 = 0;
    trans_43 = 0;
    trans_45 = 0;
    trans_54 = 0;
    
    for i=1:length(seq_idx)
        if i == length(seq_idx)
            break
        else    
            if seq_idx(i) == 1 && seq_idx(i+1) == 2
                trans_12 = trans_12 + 1;
            elseif seq_idx(i) == 2 && seq_idx(i+1) == 1
                trans_21 = trans_21 + 1;
            elseif seq_idx(i) == 2 && seq_idx(i+1) == 3
                trans_23 = trans_23 + 1;
            elseif seq_idx(i) == 3 && seq_idx(i+1) == 2
                trans_32 = trans_32 + 1;
            elseif seq_idx(i) == 3 && seq_idx(i+1) == 4
                trans_34 = trans_34 + 1;
            elseif seq_idx(i) == 4 && seq_idx(i+1) == 3
                trans_43 = trans_43 + 1;
            elseif seq_idx(i) == 4 && seq_idx(i+1) == 5
                trans_45 = trans_45 + 1;
            elseif seq_idx(i) == 5 && seq_idx(i+1) == 4
                trans_54 = trans_54 + 1; 
            end
        end
    end
    
    transitions = [trans_12, trans_23, trans_34, trans_45, trans_21, trans_32, trans_43, trans_54];
    
    Transition_prob = [];
    
    for i=1:length(transitions)
        Transition_prob = [Transition_prob; (i/(seq_num-1))*100];
    end
    
    Transition_prob = reshape(Transition_prob,[4,2]);
    
    % DISPLAY RESULTS %
    
    State = ['DFS1';'DFS2';'DFS3';'DFS4';'DFS5'];
    Transition = ['DFS1>DFS2 / DFS2>DFS1';'DFS2>DFS3 / DFS3>DFS2';'DFS3>DFS4 / DFS4>DFS3';'DFS4>DFS5 / DFS5>DFS4'];
    
    FC_dynamics1 = table(State,Fraction_times,Dwell_times);
    FC_dynamics2 = table(Transition,Transition_prob);
    
end    















    

