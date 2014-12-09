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
