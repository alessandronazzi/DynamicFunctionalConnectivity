
function [Fraction_times,Dwell_times,Transition_prob,table_seq] = FC_dyn(idx,K)
    disp('Computing clustering metrics');
    % FRACTION TIMES %
    
    Fraction_times = histcounts(idx)'/length(idx)*100;
    

    % DWELL TIMES %
   

    didx = [0 find(diff(idx))'];
    DFS_sequence=[diff(didx) length(idx)-didx(end)];
    seq_idx=idx(didx+1);
    Dwell_times = zeros(K,1);

    for i=1:K
        cc=find(seq_idx==i);
        if isempty(cc)
            Dwell_times(i)=0;
        else    
            Dwell_times(i)=mean(DFS_sequence(cc));
        end    
    end


    % TRANSITION PROBABILITY %
    
    transitions = [zeros(1,K*K)];
    seq_num = 1;
    
    for i=1:K
        for j=1:K
             if i == j
                 seq_num = seq_num + 1;
                 continue
             else    
                 for k=1:length(seq_idx)-1
                     if seq_idx(k) == i && seq_idx(k+1) == j
                        transitions(seq_num) = transitions(seq_num) + 1;
                     end   
                 end
                 seq_num = seq_num + 1;
             end    
        end
    end    
    
    
    Transition_prob = zeros(K,K);
        
    for i=1:length(transitions)
	ii=floor((i-1)/K)+1;
	jj=mod((i-1),5)+1;
        Transition_prob(ii,jj)=transitions(i)/(length(seq_idx)-1)*100;
    end
   

    Transition_prob =Transition_prob';
    
    % DISPLAY RESULTS %
    
%     State = {};
%     
%     for i=1:K
%         State{i}=['DFS',num2str(i)];
%     end    
%     
%     State = State';
%     
%     FC_dynamics1 = table(State,Fraction_times,Dwell_times);
%     FC_dynamics2 = table(Transition_prob,'RowNames',State);
end    






    

