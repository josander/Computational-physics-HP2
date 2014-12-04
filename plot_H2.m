%% Homeproblem 2b

clc
clear all

% Import data
distNuc = importdata('distances.data');
%% Plot a histogram of the distances to the nucleus

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
nbins = 50;

figure(1);
clf
[y x] = hist(distNuc, nbins);
bar(x, y/trapz(x,y))
xlabel('Distance to the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
ylabel('PDF [1/$a_0$]','Interpreter','latex', 'fontsize', 12);    
title('Probability density function for distance to nucleus','Interpreter','latex', 'fontsize', 14);
hold on
% Approximated function
r = linspace(0,4.5,100);
f = @(r) 2^5.*r.^2 .*exp(-4 .* r);
% Plot the function
plot(r,f(r), 'r', 'LineWidth', 1);
plotTickLatex2D

l = legend('Data from MC-simulation','PDF$(r) = 2^5r^2 e^{-4r}$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','distHist.eps')
%%
xlabel('Distance to the nucleus', 'fontsize', 12);
title('Approximated function', 'fontsize', 12);

%% Plot energies

% Import energy data
energy = importdata('energy.data');

% Plot data
figure(2);
clf
plot(energy(:,1), 'b');
hold on
plot(energy(:,2), 'r');
xlabel('Datapoints', 'fontsize', 12);
ylabel('Energy', 'fontsize', 12);

meanEnergy = mean(energy(:,2))


