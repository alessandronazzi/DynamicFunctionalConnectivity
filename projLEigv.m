function [upt_vmat] = projLEigv(timeseries_all,width_w,n_sbj)
%function [mat_vect] = projLEigv(timeseries_all,width_w,n_sbj)
    disp('compute pricipal eigenvectors')
    ex_check = exist('n_sbj','var');
    if ex_check == 0
        n_sbj = size(timeseries_all,3);
    end

    n_timepoints = size(timeseries_all,2);
    n_regions = size(timeseries_all,1);
    
    n_window = n_timepoints/2 - width_w + 1;

    % SLIDING-WINDOWS %

    % TO SAVE MATRICES RECONSTRUCTED FROM EIGV %
    eigv_mat = zeros(n_window*2,n_sbj,n_regions,n_regions); 
    
    % TO SAVE PRINCIPAL EIGENVECTORS %
%     mat_vect = cell(1,n_sbj); 
%     for l=1:n_sbj
%         mat_vect{l} = zeros(n_regions,n_window);
%     end            
    
    for c=1:n_sbj
        SBJn = timeseries_all(:,:,c);
      
        Mcorrs=[];
        disp(sprintf('collect sliding windows subject %d',c));   

    	top_w = 1; 
    	bot_w = width_w;
        for i=1:n_window
                timewindow_corr = corr(SBJn(:,top_w:bot_w)');
                Mcorrs = [Mcorrs; {timewindow_corr}];
                top_w = top_w+1;
                bot_w = bot_w+1;
        end
       
    	top_w = n_timepoints/2+1; 
    	bot_w = n_timepoints/2+width_w;

        for i=1:n_window
                timewindow_corr = corr(SBJn(:,top_w:bot_w)');
                Mcorrs = [Mcorrs; {timewindow_corr}];
                top_w = top_w+1;
                bot_w = bot_w+1;
        end

        
        % EIGENDECOMPOSITION AND CORRESPONDING MATRICES %
       
        disp(sprintf('compute eigendecomposition subject %d',c));   
   
        for i=1:2*n_window 
           if (mod(i,floor(2*n_window/10))==0)
               disp(i/(2*n_window)*100);
           end
           [autvet, ~, ~] = svd(Mcorrs{i});
           eigv_mat(i,c,:,:) = autvet(:,1)* transpose(autvet(:,1));
           %eigv_mat(:,i) = autvet(:,1); % TO SAVE EIGENVECTOR
        end
%         mat_vect{c} = eigv_mat; % TO SAVE EIGENVECTOR
%         eigv_mat = zeros(n_regions,n_window); % TO SAVE EIGENVECTOR
    end
       
    % UPPER-TRIANGULAR VECTORIZATION %
        
    upt_vmat = [];
    
    for c=1:n_sbj
        for i=1:2*n_window
            mat = squeeze(eigv_mat(i,c,:,:));
            upt_v = mat(triu(true(size(mat)),1));
            upt_vmat = [upt_vmat, upt_v];
        end
    end         
end    


