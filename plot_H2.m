%% Homeproblem 2b

clc
clear all

% Import data
distNuc = importdata('distances.data');

% Plot a histogram of the distances to the nucleus
nbins = 40;

figure(1);
subplot(2,1,1);
hist(distNuc, nbins);
xlabel('Distance to the nucleus', 'fontsize', 12);
title('Generated data', 'fontsize', 12);

% Approximated function
r = linspace(0,4.5,100);
f = @(r) 2^5.*r.^2 .*exp(-4 .* r);

% Plot the function
subplot(2,1,2);
plot(r,f(r));
xlabel('Distance to the nucleus', 'fontsize', 12);
title('Approximated function', 'fontsize', 12);

%% Plot energies

energy = importdata('energy.data');

figure(2);
clf
plot(energy(:,1), 'b');
hold on
plot(energy(:,2), 'r');
xlabel('Datapoints', 'fontsize', 12);
ylabel('Energy', 'fontsize', 12);
