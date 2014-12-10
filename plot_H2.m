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

print(gcf,'-depsc2','energyCorr.eps')

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

print(gcf,'-depsc2','energyBlock.eps')

%% Average energy for different alpha
a05 = [-2.86953,-2.86924,-2.87467,-2.87154,-2.86843];
a075 = [-2.87463,-2.87396,-2.86621,-2.87183,-2.87498];
a10 = [-2.87368,-2.87518,-2.87636,-2.87718,-2.87411];
a125 = [-2.87828,-2.87700,-2.87472,-2.87687,-2.87736];
a15 = [-2.87710,-2.87744,-2.87753,-2.87668,-2.87283];
a175 = [-2.87789,-2.87610,-2.87744,-2.87390,-2.87676];
a2 = [-2.87560,-2.87368,-2.87746,-2.87938,-2.87744];
a225 = [-2.87673,-2.87416,-2.87705,-2.87582,-2.87630];
a25 = [-2.87552,-2.87467,-2.87577,-2.87561,-2.87756];
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

