%%

cd 'C:\Users\aless\OneDrive\Desktop\PNC\TIRO\Schaefer2018_200Parcels_17Networks_order'
load random_silhouette_test.mat

%% SOTTOCAMPIONAMENTO %%

Npoints = [10,100,1000];
Nrand = [10000,1000,100];

sil_prova = upts';
L=size(sil_prova,1);

m_sottcamp=zeros(3,1);
std_sottcamp=zeros(3,1);

sil = [];

for i=1:3
    
    N=Npoints(i);
    %sil=[];

    for n=1:Nrand(i)
        points=datasample(1:L,N);
        sil =[sil;silhouette(sil_prova(points,:),idx(points))];
    end

    %mean_sil(i)=mean(sil,'omitnan');
    %std_sil(i)=std(sil,'omitnan');
end
%%

m_sottcamp = [mean(sil(1:100000),'omitnan'); mean(sil(100001:200000),'omitnan');
    mean(sil(200001:300000),'omitnan')];
std_sottcamp = [std(sil(1:100000),'omitnan'); std(sil(100001:200000),'omitnan'); 
    std(sil(200001:300000),'omitnan')];

errorbar(m_sottcamp,std_sottcamp)

%save('random_silhouette_test.mat','mean_sil','std_sil');




                                