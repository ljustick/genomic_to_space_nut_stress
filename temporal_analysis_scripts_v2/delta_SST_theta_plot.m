clear
close all
%This script requires scatter_kde from FileExchange

load ANOVA_monthly_SST.mat;
Season_SST = Season;

load ANOVA_monthly.mat;
Season_theta = Season;

load ANOVA_yearly_SST.mat;
Annual_SST = Annual;

load ANOVA_yearly.mat;
Annual_theta = Annual;
%% 



figure
r=1;c=2;
subplot(r,c,1)
scatter_kde(Season_SST(:),Season_theta(:));
axis([-4 4 -1 1])
xlabel('SST')
ylabel('Theta')
[Rseason,Pseason] = corr(Season_SST(:),Season_theta(:),'Rows','complete')
Rseason^2

subplot(r,c,2)
scatter_kde(Annual_SST(:),Annual_theta(:));
axis([-2 2 -1 1])
xlabel('SST')
ylabel('Theta')
[Rannual,Pannual] = corr(Annual_SST(:),Annual_theta(:),'Rows','complete')
Rannual^2
%print('-painters','-dsvg','delta_SST_theta.svg')