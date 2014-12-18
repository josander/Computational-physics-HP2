%% Homeproblem 2b

clc
clear all

% Import data
distNuc = dlmread('distances.data');
%% Plot a histogram of the distances to the nucleus
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
nbins = 150;

figure(1);
clf
[y x] = hist(distNuc, nbins);
bar(x, y/trapz(x,y))
x = xlabel('Distance from the nucleus r [$a_0$]','Interpreter','latex', 'fontsize', 12);
y = ylabel('PDF [1/$a_0$]','Interpreter','latex', 'fontsize', 12);    
title('Probability density function for distance to nucleus','Interpreter','latex', 'fontsize', 14);
hold on
% Approximated function
r = linspace(0,4.5,100);
f = @(r) 2^5.*r.^2 .*exp(-4 .* r);
% Plot the function
plot(r,f(r), 'r', 'LineWidth', 1);
axis([0 4 0 1.2])
plotTickLatex2D
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);

l = legend('Data from MC-simulation','PDF$(r) = 2^5r^2 e^{-4r}$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','distHist.eps')


%% Load energy.data

% Import energy data
energy = dlmread('energy.data');
%% Plot energies


% Plot data
figure(2);
clf
plot(energy(:,1), 'b');
hold on
plot(energy(:,2), 'r', 'LineWidth', 1.5);
x = xlabel('Iteration','Interpreter','latex', 'fontsize', 12);
y = ylabel('Energy [a.u]', 'Interpreter','latex', 'fontsize', 12);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);


meanEnergy = mean(energy(:,2))

dataS = size(energy);
block_length = 1000;

hold on

for i = 1:block_length:dataS(1) 
   plot(i+block_length - 1,mean(energy(i:i+block_length)), '. g', 'MarkerSize', 7)
  
end
axis([0 length(energy)/10 -4 -1.75])
plotTickLatex2D
l = legend('Energy','Moving energy average', 'Block averages for N = 1 000');
set(l,'Interpreter','latex')
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
%set(x, 'Units', 'Normalized', 'Position', [0.5, -0.01, 0]);
print(gcf,'-depsc2','energyAvr.eps')

%% Plot alpha

clf
plot(energy(:,3));

%%

figure(3);
plot(1:1900000,beta08, [0 1900000], [mean(beta08) mean(beta08)], 'LineWidth', 1);
x = xlabel('Iterations [ ]', 'Interpreter','latex', 'fontsize', 12);
y = ylabel('$\alpha$ [1/$a_0$]', 'Interpreter','latex', 'fontsize', 12);
axis([0 1900000 0.1 0.2])
plotTickLatex2D
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

l = legend('$\alpha$ for $\beta = 0.8$', '$\langle \alpha \rangle = 0.1482$');
set(l,'Interpreter','latex')

print(gcf,'-depsc2','alph08.eps')

%% Plot alphas for different betas
figure(6)
sizeBeta = size(beta060);


plot(1:sizeBeta(1),beta060,1:sizeBeta(1),beta070,1:sizeBeta(1),beta075,1:sizeBeta(1),beta080,1:sizeBeta(1),beta090)
axis([0 sizeBeta(1) -0.05 0.5])
plotTickLatex2D
l = legend('$\alpha$ for $\beta = 0.6$','$\alpha$ for $\beta = 0.7$','$\alpha $ for $ \beta = 0.75$','$\alpha$ for $\beta = 0.80$','$\alpha$ for $\beta = 0.9$');
set(l,'Interpreter','latex')
y = ylabel('$\alpha$ [1/$a_0$]','Interpreter','latex', 'fontsize', 12);
x = xlabel('Iteration [ ]', 'Interpreter','latex', 'fontsize', 12);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

set(y, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
print(gcf,'-depsc2','alphBet.eps')

%%
subplot(4,1,1)
plot(1:sizeBeta(1),beta060)
axis([0 sizeBeta(1) -0.05 0.5])

subplot(4,1,2)
plot(1:sizeBeta(1),beta070)
axis([0 sizeBeta(1) -0.05 0.5])

subplot(4,1,3)
plot(1:sizeBeta(1),beta075)
axis([0 sizeBeta(1) -0.05 0.5])

subplot(4,1,4)
plot(1:sizeBeta(1),beta090)
axis([0 sizeBeta(1) -0.05 0.5])

%% Corr func

% import data
data = dlmread('energy.data');
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
%subplot(1,1,1)
plot(0:numlags,corr, [0 numlags], [exp(-2) exp(-2)],'--r', i-1, corr(i),'.', 'MarkerSize', 25);
title('Auto-correlation function','Interpreter','latex','fontsize',14);
x = xlabel('Iteration lag [ ]','Interpreter','latex', 'fontsize', 12);
y = ylabel('Energy autocorrelation [ ]', 'Interpreter','latex', 'fontsize', 12);

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
l = legend('Energy autocorrelation function','$y=e^{-2}$', 'Statistical inefficiency = 11');
set(l,'Interpreter','latex')
axis([0 100 0 1.1])
plotTickLatex2D
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

print(gcf,'-depsc2','energyCorr.eps')

%% Block averaging


% import data
block = importdata('block_s.data');
blockLength = 500;

% calculate the statistical inefficiency
statistical_inefficiency = mean(block(blockLength/25:blockLength/10))

% plot
figure(5);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

plot(0:10:blockLength-10,block,'o', [0 blockLength], [statistical_inefficiency statistical_inefficiency]);
x = xlabel('Blocksize [ ]','Interpreter','latex','fontsize',12);
ylabel('Statistical inefficiency [ ]','Interpreter','latex','fontsize',12);
title('Block averaging','Interpreter','latex','fontsize',12);

set(x, 'Units', 'Normalized', 'Position', [0.5, -0.05, 0]);

plotTickLatex2D
print(gcf,'-depsc2','energyBlock.eps')

%% Average energy for different alpha
clc

a05 = [-2.86679, -2.86858, -2.86873, -2.87119, -2.86717];
a075 = [-2.87228, -2.87152, -2.87384, -2.87082, -2.87192];
a10 = [-2.87879, -2.87516, -2.87393, -2.87784, -2.87388];
a125 = [-2.87753, -2.87855, -2.87776, -2.87337, -2.87812];
a15 = [-2.87399, -2.88009, -2.87752, -2.87622, -2.87917];
a175 = [-2.87716, -2.87676, -2.87949, -2.87646, -2.87477];
a2 = [-2.87882, -2.87919, -2.87667, -2.87503, -2.87344];
a225 = [-2.87987, -2.87596, -2.87451, -2.87402, -2.87356];
a25 = [-2.87073, -2.87295, -2.87135, -2.87585, -2.87639];

mean05 = mean(a05)
mean075 = mean(a075)
mean10 = mean(a10)
mean125 = mean(a125)
mean15 = mean(a15)
mean175 = mean(a175)
mean20 = mean(a2)
mean225 = mean(a225)
mean25 = mean(a25)

var05 = var(a05)
var075 = var(a075)
var10 = var(a10)
var125 = var(a125)
var15 = var(a15)
var175 = var(a175)
var20 = var(a2)
var225 = var(a225)
var25 = var(a25)

std05 = std(a05)
%%
data = beta090;
datasq = data.^2;


norm = mean(datasq) - mean(data)^2;
cov =  xcov(data, 500000);
cov =  cov/max(cov);
index = find(cov(500001:end) < exp(-2),1);
phi_s = cov(500+index);

sigma = sqrt(norm/(900000)*index)
