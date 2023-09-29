function [timeseries_all_pr] = phase_rand(timeseries_all)
	%% X must be a matrix ntimes x nrois
  
    n_regions = size(timeseries_all,1);
    n_sbj = size(timeseries_all,3);
    n_windows = size(timeseries_all,2);
    timeseries_all_pr = zeros(n_regions,n_windows,n_sbj); 
   
    
    for s=1:n_sbj
     for sess=1:2
    	sbj = timeseries_all(:,(n_windows/2*(sess-1)+1):(n_windows/2*sess),s)';
         
        [ntimes,nroi]=size(sbj);

        Xfft = fft(sbj);

        if(rem(ntimes,2)==0)
          nph=ntimes/2-1;
          randph=2*pi*rand(nph,1);
          Xfft1=Xfft;
          Xfft1(2:nph+1,:)=Xfft(2:nph+1,:).*(exp(1i*randph)*ones(1,nroi));
          Xfft1(ntimes:-1:nph+3,:)=Xfft(ntimes:-1:nph+3,:).*(exp(-1i*randph)*ones(1,nroi));

        elseif(rem(ntimes,2)==1)
          nph=(ntimes-1)/2;
          randph=2*pi*rand(nph,1);
          Xfft1=Xfft;
          Xfft1(2:nph+1,:)=Xfft(2:nph+1,:).*(exp(1i*randph)*ones(1,nroi));
          Xfft1(ntimes:-1:nph+2)=Xfft(ntimes:-1:nph+2).*(exp(-1i*randph)*ones(1,nroi));
        end

        X1=real(ifft(Xfft1));
        X1=X1';
        timeseries_all_pr(:,(n_windows/2*(sess-1)+1):(n_windows/2*sess),s) = X1;
      end

    end
    
end
