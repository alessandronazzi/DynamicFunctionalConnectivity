function [ft,dt] = boxp_dyn(n_sbj,idx,l,width_w)

    n_window = 2*(l/2 - width_w + 1);
    ss = zeros(1, n_sbj) + n_window;
    idx_sbj = mat2cell(idx,ss,1);
    ft = [];
    dt = [];

    for i=1:length(idx_sbj)
        sbj = idx_sbj{i};
        [Fraction_times,Dwell_times,~] = FC_dyn(sbj,5);
        ft = [ft, Fraction_times];
        dt = [dt, Dwell_times];
    end

    ft = ft';
    dt = dt';
    
    figure()
    boxplot(ft)
    title('Fraction times')
    
    figure()
    boxplot(dt)
    title('Dwell times')
end    