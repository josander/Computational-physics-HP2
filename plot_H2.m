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
%%
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

% Plot data
figure(2);
clf
plot(energy(:,1), 'b');
hold on
plot(energy(:,2), 'r', 'LineWidth', 1.5);
x = xlabel('Iteration','Interpreter','latex', 'fontsize', 12);
y = ylabel('Energy [a.u]', 'Interpreter','latex', 'fontsize', 12);



meanEnergy = mean(energy(:,2))

dataS = size(energy);
block_length = 10000;

hold on

for i = 1:block_length:dataS(1) 
   plot(i+block_length - 1,mean(energy(i:i+block_length)), '. g', 'MarkerSize', 7)
  
end
axis([0 1000000 -4 -1.5])
plotTickLatex2D
l = legend('Energy','Moving energy average', 'Block averages for N = 1000');
set(l,'Interpreter','latex')
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
%set(x, 'Units', 'Normalized', 'Position', [0.5, -0.01, 0]);
print(gcf,'-depsc2','energyAvr.eps')

%% Plot alpha

figure(3);
plot(energy(:,3));
xlabel('Datapoints', 'fontsize', 12);
ylabel('Alpha', 'fontsize', 12);

%% Corr func

% import data
data = importdata('energy.data');
%%
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
plot(0:numlags,corr, [0 numlags], [exp(-2) exp(-2)],'--r', i-1, corr(i),'.', 'MarkerSize', 25);
title('Auto-correlation function','Interpreter','latex','fontsize',14);
x = xlabel('Iteration lag []','Interpreter','latex', 'fontsize', 12);
y = ylabel('Energy autocorrelation []', 'Interpreter','latex', 'fontsize', 12);



l = legend('Energy autocorrelation function','$y=e^{-2}$', 'Statistical inefficiency = 13');
set(l,'Interpreter','latex')
plotTickLatex2D
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



