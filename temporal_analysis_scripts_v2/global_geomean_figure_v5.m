% Script to generate global mean theta prime

clear
close all
addpath(genpath('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/cdt'))
addpath(genpath('/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/colorbrewer'))

load '/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/DataLL_theta.mat'
load '/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/LAT2.mat'
load '/Users/amartiny/Library/CloudStorage/GoogleDrive-amartiny@uci.edu/My Drive/NASA_NutStress/data_processed/LON3.mat'


data = DATAlowLat3;

M = geomean(data,3,'omitnan');

h = imagescn(LON2,LAT2,M);
xticks([45 90 135 180 225 270 315 360])
xticklabels({'45E','90E','135E','180E','135W','90W','45W','0'});
clim([0 4])
colorbar('eastoutside')


colors = flip(brewermap([],"Spectral"));

colormap(colors)

%print('-painters','-dsvg','global_geomean_v5.svg')

