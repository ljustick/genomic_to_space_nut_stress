% Script to generate global mean theta prime
clear
close all
%Requires CDT, inferno_map, and colorbrewer from FileExchange
addpath(genpath('cdt'))
addpath(genpath('colorbrewer'))

load inferno_cmap.mat;
load time.mat
load DataLL_theta.mat;

t = timenum;

data = [DATAlowLat3(:,201:360,:) DATAlowLat3(:,1:200,:)];

LAT = LAT(50:131);
LON2 = 20.5:379.5;
h = imagescn(LON2,LAT,mean(data,3,'omitnan'));
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
title('Annual mean NutStress')
clim([0 5])
colorbar('northoutside')


colors = flip(brewermap([],"Spectral"));
colormap(colors)

