function [DFSs] = Dfss(centroids,regions,atlas)

    disp('Visualize Dynamic Functional States')
    load('C:\Users\aless\Desktop\PNC\TIRO\Schaefer2018_200Parcels_17Networks_order\NET_info_reduction2.mat')

    DFSs = [];
    communities_gl = {'VIS';'SMN';'AUD';'CON';'VAN';'DAN';'FPN';'DMN';'none';'SUB'}; %NET_new.Names2;
    communities_sch = {'VIS';'SMN';'DAN';'VAN';'LIM';'CON';'DMN';'SUB'};
    sub_comm_rev = {'CER';'TAL';'CAU';'PUT';'PAL';'BST';'HIP';'AMY';'ACC';'DIE';'CER';'TAL'; ...
        'CAU';'PUT';'PAL';'HIP';'AMY';'ACC';'DIE'};
    sub_comm_tian = {'HIP-hm1';'HIP-hm2';'THA-VAip';'THA-VAia';'HIP-hl';'HIP-b';'HIP-t';'THA-VPm';'THA-VPl';'THA-VAs'; ...
        'THA-DAm';'THA-DAl';'PUT-VA';'PUT-DA';'PUT-VP';'PUT-DP';'CAU-VA';'CAU-DA';'CAU-b';'CAU-t';'lAMY';'mAMY';
        'THA-DP';'NAc-s';'NAc-c';'pGP';'aGP';'HIP-hm1';'HIP-hm2';'THA-VAip';'THA-VAia';'HIP-hl';'HIP-b';'HIP-t';'THA-VPm';'THA-VPl';'THA-VAs'; ...
        'THA-DAm';'THA-DAl';'PUT-VA';'PUT-DA';'PUT-VP';'PUT-DP';'CAU-VA';'CAU-DA';'CAU-b';'CAU-t';'lAMY';'mAMY';
        'THA-DP';'NAc-s';'NAc-c';'pGP';'aGP'};
    sub_comm_tian_s1 = {'HIP';'AMY';'pTHA';'aTHA';'NAc';'GP';'PUT';'CAU';'HIP';'AMY';'pTHA';'aTHA';'NAc';'GP';'PUT';'CAU'};
    reshaff_inv = [8 18 6 16 14 5 12 3 1 10 9 19 7 17 15 13 4 2 11];

    net_gl = NET_new.net_reduced;
    net_sch = [17 31 46 58 63 76 100];
    net_tian = zeros(125,1);
    net_tian_s1 = zeros(87,1);
    net_tian(1:71)=net_gl(1:71);
    net_tian_s1(1:71)=net_gl(1:71);
    net_tian(72:end)=10;
    net_tian_s1(72:end)=10;
    ticks_gl = NET_new.Tick_reduced;
    ticks_sch = [8.5,24,38.5,52,60.5,69.5,88,109.5];

    for i=1:size(centroids,1)
        C = zeros(regions);
        c = centroids(i,:);
        C(triu(true(regions),1)) = c;
        C=C+C';
        DFSs = [DFSs, {C}];
    end

    t = tiledlayout(5*2+1,5*size(centroids,1));
    
    ex_check = exist('atlas','var');

    if ex_check == 1  && strcmp(atlas,'sch')

         b=[1:9 51:58 10:15 59:66 16:23 67:73 24:30 74:78 31:33 79 80 34:37 81:89 38:50 90:100];
       
        for i=1:length(DFSs)
            name = "DFS"+i;
            dfss = DFSs{i};
            dfss(1:100,:) = dfss(b,:);
            dfss(:,1:100) = dfss(:,b);
            dfss(101:119,:) = dfss(100+reshaff_inv,:);
            dfss(:,101:119) = dfss(:,100+reshaff_inv);
    
            cm = [0 0 1; 1 1 1; 1 0 0];
            cmi = interp1([-50; 0; 50], cm, (-50:50));
    
            nexttile([5 5])
            imagesc(dfss, [-0.01 0.01])
            pbaspect([1 1 1])
            title(name)
            xticks(ticks_sch)
            yticks(ticks_sch)
            xticklabels(communities_sch)
            set(gca,'XTickLabelRotation',45,'FontSize',12)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_sch)
            else
                yticklabels([])
            end    

            for n=1:length(net_sch)

                    xline(net_sch(n)+0.5,'LineWidth',2)
                    yline(net_sch(n)+0.5,'LineWidth',2)
            end
        end
    
        %saveas(gcf,'/home/allegra/hcp_padova/results/graphs/'+name+'_'+analysistype+'.png')
    
        for i=1:length(DFSs)  
            nexttile([5 5])
            dfss = DFSs{i};
            dfss(101:119,:) = dfss(100+reshaff_inv,:);
            dfss(:,101:119) = dfss(:,100+reshaff_inv);
            imagesc(dfss(:,101:end), [-0.01 0.01])
            pbaspect([1 1 1])
            xticks(1:19)
            yticks(ticks_sch)
            xticklabels(sub_comm_rev)
            set(gca,'FontSize',10)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_sch)
            else
                yticklabels([])
            end    

            for n=1:length(net_sch) 
                    yline(net_sch(n)+0.5,'LineWidth',2)
            end
        end 
    
        for i=1:length(DFSs)
            leg = nexttile([1 5]);
            X = [0, 47, 90];
            Y = [0; 1];
            mymap = [0 0 0; 1 1 1; 0 0 0];
            C = [1 2 3; 3 3 3];
            pcolor(X,Y,C);
            pbaspect([18 1 1])
            colormap(leg,mymap);
            xticks([22 69]);
            yticklabels([])
            xticklabels({'LEFT','RIGHT'});
        end
    
        t.TileSpacing = 'compact';
        t.Padding = "tight";

    elseif ex_check==1 && strcmp(atlas,'tians4')

        for i=1:length(DFSs)
            name = "DFS"+i;
            dfss = DFSs{i};
            % dfss(72:90,:) = dfss(71+reshaff_inv,:);
            % dfss(:,72:90) = dfss(:,71+reshaff_inv);
    
            cm = [0 0 1; 1 1 1; 1 0 0];
            cmi = interp1([-50; 0; 50], cm, (-50:50));

            ticks_gl(10)=98.5;
    
            nexttile([5 5])
            imagesc(dfss, [-0.005 0.005])
            pbaspect([1 1 1])
            title(name)
            xticks(ticks_gl)
            yticks(ticks_gl)
            xticklabels(communities_gl)
            set(gca,'XTickLabelRotation',45,'FontSize',12)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_gl)
            else
                yticklabels([])
            end    
    
            for n=1:length(net_tian)
                if net_tian(n) == net_tian(end)
                    continue
                else    

                    if net_tian(n) == net_tian(n+1)
                        continue
                    else
                        xline(n+0.5,'LineWidth',2)
                        yline(n+0.5,'LineWidth',2)
                    end
                end
            end
        end
        
        for i=1:length(DFSs)  
            nexttile([5 5])
            dfss = DFSs{i};
            % dfss(101:119,:) = dfss(100+reshaff_inv,:);
            % dfss(:,101:119) = dfss(:,100+reshaff_inv);
            imagesc(dfss(:,72:end), [-0.005 0.005])
            %axis image
            pbaspect([1 1 1])
            xticks(1:54)
            yticks(ticks_gl)
            xticklabels(sub_comm_tian)
            set(gca,'FontSize',5)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_gl)
            else
                yticklabels([])
            end    

            for n=1:length(net_tian)
                if net_tian(n)==net_tian(end)
                    continue
                else    
                       
                    if net_tian(n)==net_tian(n+1)
                        continue
                    else    
                        yline(n+0.5,'LineWidth',2)
                    end
                end    
            end
        end 
    
        for i=1:length(DFSs)
            leg = nexttile([1 5]);
            X = [0, 63, 125];
            Y = [0; 1];
            mymap = [0 0 0; 1 1 1; 0 0 0];
            C = [1 2 3; 3 3 3];
            pcolor(X,Y,C);
            %axis image
            pbaspect([18 1 1])
            colormap(leg,mymap);
            xticks([31.5 94]);
            yticklabels([])
            xticklabels({'RIGHT','LEFT'});
        end
    
        t.TileSpacing = 'compact';
        t.Padding = "tight";

    elseif ex_check==1 && strcmp(atlas,'tians1')

        for i=1:length(DFSs)
            name = "DFS"+i;
            dfss = DFSs{i};
            % dfss(72:90,:) = dfss(71+reshaff_inv,:);
            % dfss(:,72:90) = dfss(:,71+reshaff_inv);
    
            cm = [0 0 1; 1 1 1; 1 0 0];
            cmi = interp1([-50; 0; 50], cm, (-50:50));

            ticks_gl(10)=79;
    
            nexttile([5 5])
            imagesc(dfss, [-0.015 0.015])
            %axis image
            pbaspect([1 1 1])
            title(name)
            xticks(ticks_gl)
            yticks(ticks_gl)
            xticklabels(communities_gl)
            set(gca,'XTickLabelRotation',45,'FontSize',12)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_gl)
            else
                yticklabels([])
            end    
    
            for n=1:length(net_tian_s1)
                if net_tian_s1(n) == net_tian_s1(end)
                    continue
                else    

                    if net_tian_s1(n) == net_tian_s1(n+1)
                        continue
                    else
                        xline(n+0.5,'LineWidth',2)
                        yline(n+0.5,'LineWidth',2)
                    end
                end
            end
        end
        
        for i=1:length(DFSs)  
            nexttile([5 5])
            dfss = DFSs{i};
            % dfss(101:119,:) = dfss(100+reshaff_inv,:);
            % dfss(:,101:119) = dfss(:,100+reshaff_inv);
            imagesc(dfss(:,72:end), [-0.015 0.015])
            %axis image
            pbaspect([1 1 1])
            xticks(1:54)
            yticks(ticks_gl)
            xticklabels(sub_comm_tian_s1)
            set(gca,'FontSize',10)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_gl)
            else
                yticklabels([])
            end    

            for n=1:length(net_tian_s1)
                if net_tian_s1(n)==net_tian_s1(end)
                    continue
                else    
                       
                    if net_tian_s1(n)==net_tian_s1(n+1)
                        continue
                    else    
                        yline(n+0.5,'LineWidth',2)
                    end
                end    
            end
        end 
    
        for i=1:length(DFSs)
            leg = nexttile([1 5]);
            X = [0, 44, 87];
            Y = [0; 1];
            mymap = [0 0 0; 1 1 1; 0 0 0];
            C = [1 2 3; 3 3 3];
            pcolor(X,Y,C);
            %axis image
            pbaspect([18 1 1])
            colormap(leg,mymap);
            xticks([22 66]);
            yticklabels([])
            xticklabels({'RIGHT','LEFT'});
        end
    
        t.TileSpacing = 'compact';
        t.Padding = "tight";
    else

        for i=1:length(DFSs)
            name = "DFS"+i;
            dfss = DFSs{i};
            dfss(72:90,:) = dfss(71+reshaff_inv,:);
            dfss(:,72:90) = dfss(:,71+reshaff_inv);
    
            cm = [0 0 1; 1 1 1; 1 0 0];
            cmi = interp1([-50; 0; 50], cm, (-50:50));
    
            nexttile([5 5])
            imagesc(dfss, [-0.015 0.015])
            %axis image
            pbaspect([1 1 1])
            title(name)
            xticks(ticks_gl)
            yticks(ticks_gl)
            xticklabels(communities_gl)
            set(gca,'XTickLabelRotation',45,'FontSize',12)
            colormap(cmi)

            if i == 1
                yticklabels(communities_gl)
            else
                yticklabels([])
            end    
            
            for n=1:length(net_gl)
                if net_gl(n) == net_gl(end)
                    continue
                else    

                    if net_gl(n) == net_gl(n+1)
                        continue
                    else
                        xline(n+0.5,'LineWidth',2)
                        yline(n+0.5,'LineWidth',2)
                    end
                end
            end
        end
    
        %saveas(gcf,'/home/allegra/hcp_padova/results/graphs/'+name+'_'+analysistype+'.png')
    
        for i=1:length(DFSs)  
            nexttile([5 5])
            dfss = DFSs{i};
            dfss(72:90,:) = dfss(71+reshaff_inv,:);
            dfss(:,72:90) = dfss(:,71+reshaff_inv);
            imagesc(dfss(:,72:end), [-0.015 0.015])
            %axis image
            pbaspect([1 1 1])
            xticks(1:19)
            yticks(ticks_gl)
            xticklabels(sub_comm_rev)
            set(gca,'FontSize',12)
            colormap(cmi)
    
            if i == 1
                yticklabels(communities_gl)
            else
                yticklabels([])
            end    
    
            for n=1:length(net_gl)
                if n == length(net_gl)
                    continue
                else
                    if net_gl(n) == net_gl(n+1)
                        continue
                    else
                        yline(n+0.5,'LineWidth',2)
                    end
                end
            end
        end 
    
        for i=1:length(DFSs)
            leg = nexttile([1 5]);
            X = [0, 47, 90];
            Y = [0; 1];
            mymap = [0 0 0; 1 1 1; 0 0 0];
            C = [1 2 3; 3 3 3];
            pcolor(X,Y,C);
            %axis image
            pbaspect([18 1 1])
            colormap(leg,mymap);
            xticks([22 69]);
            yticklabels([])
            xticklabels({'LEFT','RIGHT'});
        end
    
        t.TileSpacing = 'compact';
        t.Padding = "tight";
    
        %saveas(gcf,'/home/allegra/hcp_padova/results/graphs/'+name+'_'+'sub'+'_'+analysistype+'.png')
    end
end

