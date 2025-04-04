clear
close all
%Requires CDT and infermo_map.mat from FileExchange
addpath(genpath('cdt'))

load ANOVA_output.mat;
load DataLL_theta.mat
load time.mat
t = timenum;

LAT = LAT(50:131);
LON2 = 20.5:379.5;
data = [DATAlowLat3(:,201:360,:) DATAlowLat3(:,1:200,:)];

%This section is about detrending data
data2 = detrend3(data,t);
[eof_maps,pc,expv] = eof(data2,8);


%Plot variance explained by each EOF PC
bar(expv)
title('Variance for each EOF')
Sexpv8 = sum(expv(5:8))
Sexpv4 = sum(expv(1:4))
Sexpv1_8 = sum(expv(1:8))

%Figure with loadings
r=4;c=2;
figure
for i=1:8
    subplot(r,c,i)
    imagescn(LON2,LAT,eof_maps(:,:,i));
    colorbar
    cmocean('balance')
    clim([-0.05 0.05])
    xticks([45 90 135 180 225 270 315 360])
    xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
end

%print('-painters','-dsvg','EOF_loadings2.svg')

%Temporal evolution of each EOF
figure
r=8;c=4;
range = [-70 -60 -50 -60 -40 -40 -30 -30; 70 60 50 60 40 40 30 30];

for i=1:8
    j = 4*i-3;
    k = 4*i-1;
    ax = subplot(r,c,[j,k]);
    plot(time,pc(i,:),'-')
    axis([2002 2022 range(1,i) range(2,i)])
    grid on
    grid minor
    ax.MinorGridLineStyle = '-';
    m=4*i;
    subplot(r,c,m)
    plotpsd(pc(i,:),45)
    xlim([0 3])
    box on   
end

save('EOF_theta.mat','eof_maps','pc',"expv");