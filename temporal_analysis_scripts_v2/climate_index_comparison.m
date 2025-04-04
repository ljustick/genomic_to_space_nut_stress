clear
close all
addpath(genpath('cdt'))

load inferno_cmap.mat;
load DataLL_theta.mat;
load time.mat;
load EOF_theta.mat;

load ONI.mat;
time_oni = ONI(:,2)+ONI(:,3)/12;%ONI data was in months

load PDO.mat;
time_pdo = PDO(:,2)+PDO(:,3)/12;%ONI data was in months

figure
r=6;c=1;

YY=zeros(897,4);

for i=1:4
    j=i+4;
    YY(:,i) = smooth(pc(j,:),30,'lowess');
end


names = {'PC5','PC6','PC7','PC8'};
for i = 1:4
    subplot(r,c,i)
    plot(time,YY(:,i),'-')
    axis([2002 2022 -25 25])
    title(names(i))
    grid on 
end

subplot(r,c,5)
plot(time_oni,ONI(:,1),'-')
axis([2002 2022 -3 3])
title('ONI')
grid on


subplot(r,c,6)
plot(time_pdo,smooth(PDO(:,1),10,'loess'),'-')
axis([2002 2022 -3 3])
title('PDO')
grid on


%interpolate time for each index
%MEI
vq = zeros(897,2);
vq(:,1) = interp1(time_oni,ONI(:,1),time);
vq(:,2) = interp1(time_pdo,PDO(:,1),time);

[R,P] = corr(YY,vq);
%PC5 = ONI
%PC6 = PDO

%
figure
r=1;c=2;
subplot(r,c,1)
yyaxis left
plot(time,-YY(:,1),'-')
ylim([-25 25])
ylabel('PC5')

yyaxis right
plot(time_oni,ONI(:,1),'-')
ylim([-2.8 2.8])
ylabel('ONI')
xlim([2002 2022])
grid on

subplot(r,c,2)
yyaxis left
plot(time,YY(:,2),'-')
ylim([-25 25])
ylabel('PC6')

yyaxis right
plot(time_pdo,smooth(PDO(:,1),10,'loess'),'-')
ylim([-2.8 2.5])
ylabel('PDO')
xlim([2002 2022])
grid on
