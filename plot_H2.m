%% Homeproblem 2b

clc
clear all

nbins = 40;

% Import data
distNuc = importdata('distances.data');

% Plot a histogram of the distances to the nucleus
subplot(2,1,1);
hist(distNuc, nbins);
xlabel('Distance to the nucleus');
title('Generated data');


r = linspace(0,4.5,100);
f = @(r) 2^5.*r.^2 .*exp(-4 .* r);

subplot(2,1,2);
plot(r,f(r));
xlabel('Distance to the nucleus');
title('Approximated function');