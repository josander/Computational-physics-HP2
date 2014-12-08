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
xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
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

<<<<<<< HEAD
meanEnergy = mean(energy(:,1))
%energy(end,2)
=======
meanEnergy = mean(energy(:,2))

%% Plot alpha

figure(3);
plot(energy(:,3));
xlabel('Datapoints', 'fontsize', 12);
ylabel('Alpha', 'fontsize', 12);

%% Corr func

% import data
data = importdata('energy.data');
numlags = 200;
corr = autocorr(data(:,1), numlags);

% find statistical inefficiency
i = 1;
while corr(i) >= exp(-2)
   i = i + 1;
end

% Since no 0 index
statistical_inefficiency = i - 1

% plot
figure(4);
subplot(1,1,1)
plot(0:numlags,corr, [0 numlags], [exp(-2) exp(-2)], i, corr(i),'x');
title('Auto-correlation function','fontsize',12);
xlabel('Lags','fontsize',12);

%% Block averaging


% import data
block = importdata('block_s.data');
blockLength = 500;

% calculate the statistical inefficiency
statistical_inefficiency = mean(block(blockLength/20:blockLength/10))

% plot
figure(5);
plot(0:10:blockLength-10,block,'o', [0 blockLength], [statistical_inefficiency statistical_inefficiency]);
xlabel('Blocksize','fontsize',12);
ylabel('Statistical inefficiency','fontsize',12);
title('Block averaging','fontsize',12);

>>>>>>> 3003507b00ddb7ba6d6d9ffb3b7a373fa915430e

