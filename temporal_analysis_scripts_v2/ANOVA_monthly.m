clear
close all
%Requires CDT from FileExchange
addpath(genpath('cdt'))

load ANOVA_output.mat;
load DataLL_theta.mat

LAT = LAT(50:131);
LON2 = 20.5:379.5;

A = zeros(length(ind),12);
A(:,:) = NaN;
N = size(Coefs,2);
M = N-11;

A(ind,:) = Coefs(:,M:N);
Season = reshape(A,[82,360,12]);
Season = [Season(:,201:360,:) Season(:,1:200,:)];

Names = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
figure
r=4;c=3;
for i=1:12
    subplot(r,c,i);
    imagescn(LON2,LAT,Season(:,:,i));
    %colorbar
    clim([-2 2])
    cmocean('balance')
    xticks([45 90 135 180 225 270 315 360])
    xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
    title(Names(i))
end

%print('-painters','-dsvg','theta_monthly_variance.svg')
save("ANOVA_monthly.mat","Season")