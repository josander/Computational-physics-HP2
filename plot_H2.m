%% Homeproblem 2b

clc
clear all

nbins = 25;

% Import data
distNuc = importdata('distances.data');

% Plot a histogram of the distances to the nucleus
hist(distNuc, nbins);
xlabel('Distance to the nucleus');