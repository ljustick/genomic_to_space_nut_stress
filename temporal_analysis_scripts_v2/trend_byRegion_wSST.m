clear
close all
addpath(genpath('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/cdt'))

T_WOD = readtable("ByRegion/Dataset S1_WOD.csv");
T_WOD.lon_degE = T_WOD.lon_degE+0.5;
T_WOD.lat_degN = T_WOD.lat_degN-0.5;

T_WOD = T_WOD(abs(T_WOD.lat_degN)< 40.6,:);

load('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/DataLL_theta.mat');
load('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/time.mat');

theta = DATAlowLat3;
Tsq = 365.25*trend(theta,timenum,3);
T = Tsq(:);

load('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/DataLL_SST.mat');
T_SST_sq = 365.25*trend(DATAlowLat3,timenum,3);
T_SST = T_SST_sq(:);

load('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/LON2.mat');
LON = LON2;
load('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/LAT.mat');
LAT = LAT(50:131);

%%
A = [T_WOD.lat_degN T_WOD.lon_degE];
MLON = repmat(LON',82,1);
MLAT = repmat(LAT',360,1);
MLAT = MLAT';
Grid_theta = [MLAT(:) MLON(:)];

[C,ia,ib] = intersect(A,Grid_theta,'rows');

T_WOD.theta_trend(ia) = T(ib);
T_WOD.SST_trend(ia) = T_SST(ib);

T_WOD = T_WOD(~isnan(T_WOD.theta_trend),:);

%%
r = 4;c=6;
edges = -0.03:0.001:0.03;

subplot(r,c,1:2)
histogram(T_WOD.theta_trend,edges,'EdgeColor','none')
title('all')
xline(0)

m = bootstrp(10000,@mean,T_WOD.theta_trend);
[fi,xi] = ksdensity(m);
subplot(r,c,3)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

subplot(r,c,4:5)
histogram(T_WOD.theta_trend,edges,'EdgeColor','none')
title('all')
xline(0)

m = bootstrp(10000,@mean,T_WOD.theta_trend);
[fi,xi] = ksdensity(m);
subplot(r,c,6)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)


subplot(r,c,7:8)
ind = find(T_WOD.Zno3_m>50);
histogram(T_WOD.theta_trend(ind),edges,'EdgeColor','none')
title('All >50 m')
xline(0)

m = bootstrp(10000,@mean,T_WOD.theta_trend(ind));
[fi,xi] = ksdensity(m);
subplot(r,c,9)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

subplot(r,c,10:11)
ind = find(T_WOD.Zno3_m<50);
histogram(T_WOD.theta_trend(ind),edges,'EdgeColor','none')
title('All <50 m')
xline(0)

m = bootstrp(10000,@mean,T_WOD.theta_trend(ind));
[fi,xi] = ksdensity(m);
subplot(r,c,12)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

subplot(r,c,13:14)
ind = find(T_WOD.Zno3_m>50 & T_WOD.lat_degN>0);
histogram(T_WOD.theta_trend(ind),edges,'EdgeColor','none')
title('NH >50 m')
xline(0)

m = bootstrp(10000,@mean,T_WOD.theta_trend(ind));
[fi,xi] = ksdensity(m);
subplot(r,c,15)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

subplot(r,c,16:17)
ind = find(T_WOD.Zno3_m<50 & T_WOD.lat_degN>0);
histogram(T_WOD.theta_trend(ind),edges,'EdgeColor','none')
title('NH <50 m')
xline(0)

m = bootstrp(10000,@mean,T_WOD.theta_trend(ind));
[fi,xi] = ksdensity(m);
subplot(r,c,18)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

subplot(r,c,19:20)
ind = find(T_WOD.Zno3_m>50 & T_WOD.lat_degN<0);
histogram(T_WOD.theta_trend(ind),edges,'EdgeColor','none')
title('SH >50 m')
xline(0)
m = bootstrp(10000,@mean,T_WOD.theta_trend(ind));
[fi,xi] = ksdensity(m);
subplot(r,c,21)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

subplot(r,c,22:23)
ind = find(T_WOD.Zno3_m<50 & T_WOD.lat_degN<0);
histogram(T_WOD.theta_trend(ind),edges,'EdgeColor','none')
title('SH <50 m')
xline(0)
m = bootstrp(10000,@mean,T_WOD.theta_trend(ind));
[fi,xi] = ksdensity(m);
subplot(r,c,24)
plot(xi,fi)
xlim([-0.01 0.01])
xline(0)

%print('-painters','-dsvg','theta_trend_region_v3.svg')
%%


%All
i = 1;
result(1,i)     = mean(T_WOD.theta_trend);
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend);

result2(1,i)     = mean(T_WOD.SST_trend);
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend);

%All nutricline
i = i+1;
ind = find(T_WOD.Zno3_m>50);
result(1,i)     = mean(T_WOD.theta_trend(ind));
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend(ind));

result2(1,i)     = mean(T_WOD.SST_trend(ind));
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend(ind));


i = i+1;
ind = find(T_WOD.Zno3_m<50);
result(1,i)     = mean(T_WOD.theta_trend(ind));
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend(ind));

result2(1,i)     = mean(T_WOD.SST_trend(ind));
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend(ind));

%NH
i = i+1;
ind = find(T_WOD.Zno3_m>50 & T_WOD.lat_degN>0);
result(1,i)     = mean(T_WOD.theta_trend(ind));
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend(ind));

result2(1,i)     = mean(T_WOD.SST_trend(ind));
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend(ind));

i = i+1;
ind = find(T_WOD.Zno3_m<50 & T_WOD.lat_degN>0);
result(1,i)     = mean(T_WOD.theta_trend(ind));
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend(ind));

result2(1,i)     = mean(T_WOD.SST_trend(ind));
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend(ind));

%SH
i = i+1;
ind = find(T_WOD.Zno3_m>50 & T_WOD.lat_degN<0);
result(1,i)     = mean(T_WOD.theta_trend(ind));
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend(ind));

result2(1,i)     = mean(T_WOD.SST_trend(ind));
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend(ind));

i = i+1;
ind = find(T_WOD.Zno3_m<50 & T_WOD.lat_degN<0);
result(1,i)     = mean(T_WOD.theta_trend(ind));
result(2:3,i)   = bootci(10000,@mean,T_WOD.theta_trend(ind));

result2(1,i)     = mean(T_WOD.SST_trend(ind));
result2(2:3,i)   = bootci(10000,@mean,T_WOD.SST_trend(ind));

result(4,:)  = result(1,:)-result(2,:);
result(5,:)  = result(3,:)-result(1,:);

result2(4,:)  = result2(1,:)-result2(2,:);
result2(5,:)  = result2(3,:)-result2(1,:);

%%
figure
result = 10*result;
errorbar(1:7,result2(1,:)',result2(4,:)',result2(5,:)','s',"MarkerSize",10,...
    "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],"CapSize",0);

hold on

errorbar(1:7,result(1,:)',result(4,:)',result(5,:)','s',"MarkerSize",10,...
    "MarkerEdgeColor","red","MarkerFaceColor","red","CapSize",0);

ylim([-0.06 0.06])
xlim([0 8])
yline(0)
camroll(-90)

%print('-painters','-dsvg','theta_trend_region_errorbar_SST_v3.svg')