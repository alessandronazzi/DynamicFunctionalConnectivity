n_sbj = 80; %size(results.timeseries_all,3)
    n_window = 773;
    for i=1:80
        mat_vect{i}=mat_vect{i}';
    end    
    net = NET_new.net_reduced;
    names = NET_new.Names2;
    names_n = cellstr(names);
    names_n{10} = 'Null';
    %upts = results.upts(1:n_window*6,:)
    
%     net = zeros(119,1);
%     a = [17 31 46 58 63 76 100 119];
%     b = [1 18 32 47 59 64 77 101];
%     for i=1:length(a)
%         net(b(i):a(i))=i;
%     end
%     
%     names = {'VIS';'SMN';'DAN';'VAN';'LIM';'CON';'DMN';'SUB'};
%     names_n = {'VIS';'SMN';'DAN';'VAN';'LIM';'CON';'DMN';'null'};

    D_nets = zeros(10,n_window-1,n_sbj);
    nets = [];

    for s=1:n_sbj
        sbj = mat_vect{s};
        for n=1:10
            netidx = find(net==n);
            for w=1:n_window
                M=sbj(w,netidx);
                
                nets = [nets, {M}];
            end
            for i=1:length(nets)
                if i == length(nets)
                    continue
                else
                    d_net = mean(abs(nets{i+1}-nets{i}));
                    D_nets(n,i,s) = d_net;
                end
            end
            nets = [];
        end                
    end 
    
    edges = [];
    binw = [];
    values = [];

    for r=1:size(D_nets,1)
        h = histogram(D_nets(r,:,:),'Normalization','probability');
        edges = [edges,{h.BinEdges}];
        binw = [binw,{h.BinWidth}];
        values = [values,{h.Values}];
    end 

    close all
    
    subplot(2,1,1)
    for p=1:size(D_nets,1)
        hold on
        plot(edges{p}(1:end-1)+binw{p},values{p})
    end
    
    prct_95 = prctile(reshape(D_nets,nnz(D_nets),1),95);
    xline(prct_95,'LineWidth',1,'LineStyle','--','Color','r','HandleVisibility','off')
    legend(cellstr(names));
    xlabel('|DFC(t+1)-DFC(t)|');
    ylabel('Probability');
    
    conn_jumps = zeros(size(D_nets,1)-1,n_sbj);

    for s=1:size(D_nets,3)
        D_nets_n = D_nets(:,:,s);
        for n=1:size(D_nets_n,1)-1
            for w=1:size(D_nets_n,2)
                if D_nets_n(n,w) >= prct_95 && D_nets_n(10,w) >= prct_95
                    conn_jumps(n,s) = conn_jumps(n,s)+1;
                else
                    continue
                end
            end
        end
    end
    
    cort_jumps = zeros(size(D_nets,1)-1,n_sbj);
    
     for s=1:size(D_nets,3)
        D_nets_n = D_nets(:,:,s);
        for n=1:size(D_nets_n,1)-1
            for w=1:size(D_nets_n,2)
                if D_nets_n(n,w) >= prct_95 
                    cort_jumps(n,s) = cort_jumps(n,s)+1;
                else
                    continue
                end
            end
        end
     end
    
    sub_jumps = zeros(1,n_sbj);
    
    for s=1:size(D_nets,3)
        D_nets_n = D_nets(10,:,s);
        for w=1:size(D_nets_n,2)
            if D_nets_n(w) >= prct_95 
                sub_jumps(s) = sub_jumps(s)+1;
            else
                continue
            end
        end
    end
    
    for i=1:size(cort_jumps,1)
        cort_jumps(i,:) = cort_jumps(i,:)/size(D_nets,2);
    end    
    
    for i=1:length(sub_jumps)
        sub_jumps(i) = sub_jumps(i)/size(D_nets,2);
    end         

    for i=1:size(conn_jumps,1)
        conn_jumps(i,:)=(conn_jumps(i,:)/size(D_nets,2));
    end
    
    for i=1:size(conn_jumps,1)
        conn_jumps(i,:)=conn_jumps(i,:)./cort_jumps(i,:);
    end    
        
    hold off
    subplot(2,1,2)
    hold on
    
    for g=1:size(conn_jumps,1)
        %cdfplot(conn_jumps(g,:))
        plot(sort(conn_jumps(g,:)),(1:size(D_nets,3))*1/size(D_nets,3),'LineWidth',3)
    end
    
    %cdfplot(sub_jumps)
    plot(sort(sub_jumps),(1:size(D_nets,3))*1/size(D_nets,3),'LineWidth',3,'Color','k')
    
    legend(names_n)
    xlabel('P(subcortical change|cortical change)');
    ylabel('Cumulative distribution');
