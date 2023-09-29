Npoints = 70;
    names = {'VIS';'SMN';'AUD';'CON';'VAN';'DAN';'FPN';'DMN';'none'}; 
     
    ss = zeros(1,n_sbj) + n_window;
    idx_sbj = mat2cell(idx,ss,1);

    SN = zeros(11,Npoints);
    sbj_n = mat_vect{sbj};
    idx_n = idx_sbj{sbj};

    for t=strt:strt+Npoints
        M=sbj_n(t,:);
        for n=1:9
            netidx = find(net==n);
            net_n = M(netidx);
            
            %SN(n,t)=mean(mean(M(netidx,netidx)));
            SN(n,t)=mean(net_n);
        end
        netidx = find(net==10);
        rs_netidx = [netidx(1:2);netidx(5:11);netidx(14:19);netidx(3:4);netidx(12:13)];
        SN(10,t)=mean(M(rs_netidx(1:15)));
        SN(11,t)=mean(M(rs_netidx(16:19)));
    end

    %ticks = [];
    %tiledlayout(2,1)
    %nexttile
    ax1 = subplot(2,1,1);
    plot(SN(1:9,:)');
    xlim([strt strt+Npoints])
    title('Cortical networks')
    legend(ax1,cellstr(names),'Location','eastoutside');

    for i=strt:strt+Npoints
        if i == strt+Npoints
            continue
        else
            if idx_n(i)==idx_n(i+1)
                continue
            else
                xline(i+0.5,'LineWidth',1,'LineStyle','--','HandleVisibility','off')
            end
        end
    end    

    %nexttile
    ax2 = subplot(2,1,2);
    plot(SN(10:11,:)')
    xlim([strt strt+Npoints])
    title('Sub-cortical networks')
    legend(ax2,'SC-1','SC-2','Location','eastoutside')

    for i=strt:strt+Npoints
        if i == strt+Npoints
            continue
        else
            if idx_n(i)==idx_n(i+1)
                continue
            else
                xline(i+0.5,'LineWidth',1,'LineStyle','--','HandleVisibility','off')
            end
        end
    end    
