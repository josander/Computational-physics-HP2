%% Homeproblem 2b

clc
clear all

% Import data
distNuc = importdata('distances.data');

% Plot a histogram of the distances to the nucleus
hist(distNuc);