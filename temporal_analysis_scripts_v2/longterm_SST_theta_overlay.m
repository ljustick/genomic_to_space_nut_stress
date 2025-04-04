clear
close all
addpath(genpath('cdt'))
addpath(genpath('colorbrewer'))
load time.mat;



load DataLL_SST.mat;
data_sst = [DATAlowLat3(:,201:360,:) DATAlowLat3(:,1:200,:)];

load DataLL_theta.mat;
data = [DATAlowLat3(:,201:360,:) DATAlowLat3(:,1:200,:)];

LON2 = 20.5:379.5;
LAT = LAT(50:131);
T = 365.25*trend(data,timenum,3);
T_sst = 365.25*trend(data_sst,timenum,3);


edges = -0.1:0.001:0.1;
r=1;c=2;
subplot(r,c,1)
imagescn(LON2,LAT,T)
colorbar
cmocean('balance')
title('Long-time Trend NutStress')
caxis([-0.05 0.05])
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});


subplot(r,c,2)
imagescn(LON2,LAT,T_sst)
colorbar
cmocean('balance')
title('Long-time Trend SST')
caxis([-0.1 0.1])
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});


[R,P] = corr(T(:),T_sst(:),'Rows','complete')
%% Printing function
%print('-painters','-dsvg','longterm_SST_theta.svg')


%% Overlay of the two trends
%total non-NaN -> 22468
cmap_overlay = [10 70 130; 255 205 0]/255;%Maroon, Yellow


T_sst(isnan(T)) = NaN;
ind1 = find(T_sst>0 & T>0);
ind2 = find(T_sst<0 & T>0);
ind3 = find(T_sst>0 & T<0);
ind4 = find(T_sst<0 & T<0);

Toverlay = zeros(82,360);
Toverlay(isnan(T)) = NaN;
Toverlay([ind1' ind4']) = 1;
Toverlay([ind2' ind3']) = 2;



figure
imagescn(LON2,LAT,Toverlay)
%colorbar
colormap(cmap_overlay)
title('Long-time Trend Theta vs SST')
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});

%Check for missing points
%Sum = length(ind1)+ length(ind2)+length(ind3)+length(ind4);
%print('-painters','-dsvg','longterm_overlay.svg')
%% 
R1N = length(find(Toverlay(1:41,:)==1))/length(find(~isnan(Toverlay(1:41,:))))
R1S = length(find(Toverlay(42:82,:)==1))/length(find(~isnan(Toverlay(42:82,:))))
%% Area weighted values
[lat,lon] = cdtgrid(1);
A = cdtarea(lat,lon,'km^2');
A = [A(50:131,201:360,:) A(50:131,1:200,:)];
AN = A(1:41,:);
AS = A(42:82,:);

I1 = find(Toverlay(1:41,:) ==1);
I2 = find(Toverlay(1:41,:) ==2);
A1N = AN(I1);
A2N = AN(I2);

I1 = find(Toverlay(42:82,:) ==1);
I2 = find(Toverlay(42:82,:) ==2);
A1S = AS(I1);
A2S = AS(I2);
R1NA = sum(A1N)/(sum(A1N)+sum(A2N))
R2NA = sum(A2N)/(sum(A1N)+sum(A2N))

R1SA = sum(A1S)/(sum(A1S)+sum(A2S))
R2SA = sum(A2S)/(sum(A1S)+sum(A2S))
%% Area weighted for theta and SST
%Proportions
ind_high = find(T(:)>0);
ind_low  = find(T(:)<0);
R_theta = sum(A(ind_high))/(sum(A(ind_high))+sum(A(ind_low)))

ind_high = find(T_sst(:)>0);
ind_low  = find(T_sst(:)<0);
R_sst = sum(A(ind_high))/(sum(A(ind_high))+sum(A(ind_low)))
%Means
M_theta = mean(T(~isnan(T)).*A(~isnan(T)))/mean(A(~isnan(T)));
Median_theta = median(T(~isnan(T)));
M_SST   = mean(T_sst(~isnan(T_sst)).*A(~isnan(T_sst)))/mean(A(~isnan(T_sst)));
Median_SST = median(T_sst(~isnan(T_sst)));
SST.period = M_SST*(time(end) - time(1));
