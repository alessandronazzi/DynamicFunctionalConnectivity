
function [timeseries_all] = preprocess(varargin)

    nargin=length(varargin);
    timeseries_all=varargin{1};

    atlas=varargin{2};

    if(strcmp(atlas,'GordonLaumann') || strcmp(atlas,'TianS1') || strcmp(atlas,'TianS4') || strcmp(atlas,'Buckner') || strcmp(atlas,'Morel') )
        load('/home/allegra/hcp_padova/results/GordonLaumann/NET_info_reduction2.mat');
    end
    

    glob=true;
    for i=1:nargin
        if ischar(varargin{i}) && strcmp(varargin{i},'noglob')
	   glob=false;
	end
    end


    % DETREND % 

    timeseries = [];
    w_timeseries = size(timeseries_all,1);
    l_timeseries = size(timeseries_all,2);
    if length(size(timeseries_all))>2
        n_sbj = size(timeseries_all,3);
    else
        n_sbj = 1;
    end    
    
    timeseries_det = zeros(w_timeseries,l_timeseries,n_sbj);
    
    for i=1:size(timeseries_all,3)
        sbj_n = timeseries_all(:,:,i);
        sbj_n = [detrend(sbj_n(:,1:l_timeseries/2)')',detrend(sbj_n(:,l_timeseries/2+1:l_timeseries)')'];
        timeseries_det(:,:,i) = sbj_n;
    end    

    timeseries_all = timeseries_det;


    % BAND-PASS FILTER % 
    
    bpss_low=0.009;
    bpss_high=0.08;
    TR=0.71;
    filterorder=1;
    lopass=bpss_low/(0.5/TR);
    hipass=bpss_high/(0.5/TR);
    
    timeseries_lpf = zeros(w_timeseries,l_timeseries,n_sbj);
    
    for i=1:size(timeseries_all,3) 
        sbj_n = timeseries_all(:,:,i);

        [butta, buttb] = butter(filterorder,[lopass hipass]);
        sbj_n = [filtfilt(butta,buttb,sbj_n(:,1:l_timeseries/2)')',filtfilt(butta,buttb,sbj_n(:,l_timeseries/2+1:l_timeseries)')'];
        timeseries_lpf(:,:,i) = sbj_n;

    end


    timeseries_all = timeseries_lpf;

    % GLOBAL SIGNAL REGRESSION %

    if(glob)
        timeseries_noG = zeros(w_timeseries,l_timeseries,n_sbj);

        for i=1:size(timeseries_all,3)
            sbj = timeseries_all(:,:,i);
            global_signal = mean(sbj,1);
            beta_g = sbj*pinv(global_signal);
            res_sbj = sbj - beta_g*global_signal;

            timeseries_noG(:,:,i) = res_sbj;
        end
        timeseries_all = timeseries_noG;
    end

    % DIMENSIONALITY REDUCTION: 90X90 % 

    if strcmp(atlas,'GordonLaumann')

        atlas_90 = NET_new.final_index;
        conv_areas = [];
	
	timeseries_red = zeros(size(timeseries_all,1)-253,size(timeseries_all,2),size(timeseries_all,3));

        for i=1:size(timeseries_all,3)
            sbj = timeseries_all(:,:,i);
            for r=1:length(atlas_90)
                reg = atlas_90{r};
                for c=1:length(reg)
                    pos = reg(c);
                    conv_areas = [conv_areas; sbj(pos,:)];
                end
                if size(conv_areas,1) > 1
                    new_area = mean(conv_areas);
                    timeseries_red(r,:,i) = new_area;
                else
                    timeseries_red(r,:,i) = conv_areas;
                end
                conv_areas = [];  
            end

        end  
         
	timeseries_all = timeseries_red;
    
    elseif strcmp(atlas,'TianS1') || strcmp(atlas,'TianS4') || strcmp(atlas,'Buckner') || strcmp(atlas,'Morel') 
        atlas_90 = NET_new.final_index(1:71);
   
        conv_areas = [];
	
	timeseries_red = zeros(size(timeseries_all,1)-253,size(timeseries_all,2),size(timeseries_all,3));

        for i=1:size(timeseries_all,3)
            sbj = timeseries_all(:,:,i);
            for r=1:length(atlas_90)
                reg = atlas_90{r};
                for c=1:length(reg)
                    pos = reg(c);
                    conv_areas = [conv_areas; sbj(pos,:)];
                end
                if size(conv_areas,1) > 1
                    new_area = mean(conv_areas);
                    timeseries_red(r,:,i) = new_area;
                else
                    timeseries_red(r,:,i) = conv_areas;
                end
                conv_areas = [];  
            end

	    sz=size(timeseries_all,1)-324;

	    for r=1:sz
		timeseries_red(71+r,:,i) = sbj(324+r,:);
            end
        end  
        timeseries_all = timeseries_red;
    end
end
      
