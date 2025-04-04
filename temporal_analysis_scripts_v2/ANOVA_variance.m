clear
close all
%Requires CDT and infermo_map.mat from FileExchange
addpath(genpath('cdt'))

load inferno_cmap.mat;
load ANOVA_output.mat;
load DataLL_theta.mat;

LAT = LAT(50:131);
LON2 = 20.5:379.5;

A = zeros(length(ind),1);
A(:,:) = NaN;
A(ind,:) = SSQ(:,1);
SumSSQ.annual = reshape(A,[82,360]);
SumSSQ.annual = [SumSSQ.annual(:,201:360,:) SumSSQ.annual(:,1:200,:)];



A = zeros(length(ind),1);
A(:,:) = NaN;
A(ind,:) = SSQ(:,1)./SSQ(:,4);
SumSSQ.frac_annual = reshape(A,[82,360]);
SumSSQ.frac_annual = [SumSSQ.frac_annual(:,201:360,:) SumSSQ.frac_annual(:,1:200,:)];

A = zeros(length(ind),1);
A(:,:) = NaN;
A(ind,:) = SSQ(:,2);

SumSSQ.season = reshape(A,[82,360]);
SumSSQ.season = [SumSSQ.season(:,201:360,:) SumSSQ.season(:,1:200,:)];

A = zeros(length(ind),1);
A(:,:) = NaN;
A(ind,:) = SSQ(:,2)./SSQ(:,4);
SumSSQ.frac_season = reshape(A,[82,360]);
SumSSQ.frac_season = [SumSSQ.frac_season(:,201:360,:) SumSSQ.frac_season(:,1:200,:)];

%Variance plot
figure
r=2;c=2;


title('Annual SSQ NutStress')

subplot(r,c,1);
imagescn(LON2,LAT,SumSSQ.frac_annual);
colorbar
colormap(inferno_cmap)
clim([0 0.20])
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
title('Annual SSQ NutStress')


subplot(r,c,3);
imagescn(LON2,LAT,SumSSQ.frac_season);
colorbar
clim([0 0.5])
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
title('Seasonal SSQ NutStress')

subplot(r,c,[2 4]);
histogram(SumSSQ.frac_annual,'EdgeColor','none')
hold on
histogram(SumSSQ.frac_season,'EdgeColor','none')
hold off
axis([0 0.6 0 2500])
xticks(0:0.2:0.6)
yticks([0:500:2500]);
legend({'Annual','Seasonal'})

Mf_season = median(SumSSQ.frac_season(:),'omitnan')
Mf_annual = median(SumSSQ.frac_annual(:),'omitnan')