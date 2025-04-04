clear
close all
addpath(genpath('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/Inpaint_nans'))
load '/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/raw_data/thetaPrime_Dec24.mat';
load '/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed'/BathyLL.mat
DATA = double(thetaprime2);
DATA = [DATA(:,201:360,:) DATA(:,1:200,:)];

Nmax = 4;

%Low latitude array
%LAT -40 (row 131) - 40 (row 50)
DATAlowLat = DATA(50:131,:,:);
DATAlowLat2 = reshape(DATAlowLat,29520,897);
Data2 = zeros(29520,897);
Data2(Data2==0) = NaN;
Data2(bathyLL(:),:) = DATAlowLat2(bathyLL(:),:);
DATAlowLat2 = Data2;



%45 band
%Lat2 -45 (row 45) to 45 (row 136)
%DATAlowLat = DATA(45:136,:,:);
%DATAlowLat2 = reshape(DATAlowLat,[],828);


%Remove positions with >25 percent NaNs
C = sum(isnan(DATAlowLat2),2);
DATAlowLat2(C>224,:) = NaN;

figure
subplot(1,2,1)
histogram(DATAlowLat2)




%Removing outliers
DATAlowLat2(DATAlowLat2>Nmax) = Nmax;
DATAlowLat2(DATAlowLat2<0) = 0;

%Filling out NaNs with median values prior to smoothing
C = sum(isnan(DATAlowLat2),2);
ind = find(C<897 & C>0);

temp = zeros(length(ind),897);
for j = 1:length(ind)
    temp(j,:) = inpaint_nans(DATAlowLat2(ind(j),:),1);
end


DATAlowLat2(ind,:) = temp;
C2 = sum(isnan(DATAlowLat2),2);
ind2 = find(C2<897 & C2>0);


%Smooth of data (also removes NaNs)
D = DATAlowLat2';
N = size(D,2);
M = size(D,1);
DATA_LL_smooth = zeros(size(D));

for i =1:N
    if C(i) == M
        DATA_LL_smooth(:,i) = NaN;
    else
        DATA_LL_smooth(:,i) = smooth(D(:,i),5,'lowess');
    end
end

DATA_LL_smooth_temp = DATA_LL_smooth';
DATAlowLat3 = reshape(DATA_LL_smooth_temp,[82,360,897]);%82 rows if 40 band, 92 if 45 band
DATAlowLat3(DATAlowLat3>Nmax) = Nmax;
DATAlowLat3(DATAlowLat3<0) = 0;
figure
histogram(DATAlowLat3)



%save('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/DataLL_theta.mat','DATAlowLat3')
