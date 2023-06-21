clear
close all

%Requires CDT from FileExchange
addpath(genpath('cdt'))
load ANOVA_output.mat;
load DataLL_theta.mat

LAT = LAT(50:131);
LON2 = 20.5:379.5;



N = size(Coefs,2)-13;
M = N+1;
A = zeros(length(ind),N);
A(:,:) = NaN;
A(ind,:) = Coefs(:,2:M);
Annual = reshape(A,[82,360,N]);
Annual = [Annual(:,201:360,:) Annual(:,1:200,:)];
X = 2002:2021;
figure
r=5;c=4;
for i=1:N
    subplot(r,c,i);
    
    imagescn(LON2,LAT,Annual(:,:,i));
    
    colorbar
    clim([-1 1])
    cmocean('balance')
    xticks([45 90 135 180 225 270 315 360])
    xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
    title(num2str(X(i)))
end

%print('-painters','-dsvg','theta_yearly_variance.svg')
save("ANOVA_yearly.mat","Annual")