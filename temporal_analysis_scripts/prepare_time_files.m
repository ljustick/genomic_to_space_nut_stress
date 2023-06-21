clear
close all
%Requires ddd2mmdd from FileExchange

load thPrime_1degree_032223.mat;
addpath(genpath('ddd2mmdd'))

[M,D] = ddd2mmdd(yr,yd);
time_vec = [yr M D];
timenum = datenum(yr,M,D);
time = yr + yd/365;
timefile = 'time.mat';
save(timefile,'time_vec','timenum','time');