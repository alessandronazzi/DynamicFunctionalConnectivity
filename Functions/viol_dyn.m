function [ft,dt] = viol_dyn(n_sbj,idx,l,width_w,K)

    n_window = 2*(l/2 - width_w + 1);
    ss = zeros(1, n_sbj) + n_window;
    idx_sbj = mat2cell(idx,ss,1);
    ft = [];
    dt = [];

    for i=1:length(idx_sbj)
        sbj = idx_sbj{i};
        [Fraction_times,Dwell_times,~] = FC_dyn(sbj,K);
        ft = [ft, Fraction_times];
        dt = [dt, Dwell_times];
    end

    ft = ft';
    dt = dt*2.1;
    dt = dt';
    
    ft([236 487],:)=[];
    dt([236 487],:)=[];
    
    figure()
    violin(ft,'facecolor',[0.5 0.7 1],'medc',[])
    title('Fraction times')
    xlabel('DFS')
    xticks([1 2])
    ylabel('%')
    set(gca,'FontSize',15)
    
    figure()
    violin(dt,'facecolor',[0.5 0.7 1],'medc',[])
    title('Dwell times')
    xlabel('DFS')
    xticks([1 2])
    ylabel('seconds')
    set(gca,'FontSize',15)
end    
