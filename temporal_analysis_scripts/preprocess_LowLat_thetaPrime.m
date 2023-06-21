clear
close all

%Requires Inpaint_nans.m from FileExchange

addpath(genpath('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/Inpaint_nans'))
load thPrime_1degree_032223.mat;
DATA = double(thetaprime);

%Low latitude array
%LAT -40 (row 131) to 40 (row 50)
DATAlowLat = DATA(50:131,:,:);
DATAlowLat2 = reshape(DATAlowLat,29520,897);


%Remove positions with >25 percent NaNs
C = sum(isnan(DATAlowLat2),2);
DATAlowLat2(C>224,:) = NaN;

% figure
% subplot(1,2,1)
% histogram(DATAlowLat2)




%Removing outliers
DATAlowLat2(DATAlowLat2>6) = 6;
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


D = DATAlowLat2';
N = size(D,2);
M = size(D,1);
DATA_LL_smooth = zeros(size(D));

for i =1:N
    if C(i) == M
        DATA_LL_smooth(:,i) = NaN;
    else
        DATA_LL_smooth(:,i) = smooth(D(:,i),5,'loess');
    end
end

DATA_LL_smooth_temp = DATA_LL_smooth';
DATAlowLat3 = reshape(DATA_LL_smooth_temp,[82,360,897]);%82 rows if 40 band,

save('DataLL_theta.mat',"DATAlowLat3","LAT","LON")

% subplot(1,2,2)
% histogram(DATAlowLat3)
