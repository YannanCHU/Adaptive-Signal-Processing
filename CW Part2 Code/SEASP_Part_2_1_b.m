% Implement an LMS adaptive predictor using 1000 samples
% plot the learning curve

clc;
clear;
close all;

% number of samples
N = 1000;
% step sizes
mu1 = 0.05;
mu2 = 0.01;
mus = [mu1, mu2];
% coefficients of AR model
a1 = 0.1;
a2 = 0.8;
% noise variance
noiseVar = 0.25;

%% single realisation
% get the AR(2) model
% x(n) = 0 + a1*x(n-1) + a2*x(n-2) + noise
AR_MDL = arima('Constant',0,'AR',{a1 a2},'Variance',noiseVar);
% generate the siganl with N samples
rng(0);
x1 = simulate(AR_MDL, N, 'NumPaths', 1);   % simulate AR2 process

x1_errors_sq = cell(2,1);
for i = 1:length(mus)
    [~, x1_errors, ~] = LMS_arma([], x1, mus(i), 2, 0);
    x1_errors_sq{i} = 10 * log10(x1_errors .^ 2);
end

figure();
tiledlayout(2,1,'TileSpacing','compact');   nexttile;
plot(1:N, cell2mat(x1_errors_sq.').', 'LineWidth', 2);
title("The learning curve of one realisation of x(n)", 'FontSize',16);
ylabel("Squared prediction error in dB", 'FontSize',14);
xlabel("Time index", 'FontSize',14);
legend("\mu = "+ mu1, "\mu = " + mu2, 'FontSize',12, 'Location','southeast');


%% 100 realisations and the averaged learning curve
% generate 100 realisations of x(n)
rng(0);
numOfRealisations = 100;
x100 = simulate(AR_MDL, N, 'NumPaths', numOfRealisations);   % simulate AR2 process
x100_errors_sq = cell(length(mus),numOfRealisations);         % 2 by 100 cell

for j = 1:numOfRealisations
    
    for i = 1:length(mus)
        [~, errors, ~] = LMS_arma([], x100(:,j), mus(i), 2, 0);
        x100_errors_sq{i, j} = 10 * log10(errors .^ 2);
    end
end

nexttile; plot(1:N, mean(cell2mat(x100_errors_sq(1,:)).'), ...
    1:N, mean(cell2mat(x100_errors_sq(2,:)).'),'LineWidth', 2);
title("The aveage learning curve of 100 realisations of x(n)", 'FontSize',16);
ylabel("Squared prediction error in dB", 'FontSize',14);
xlabel("Time index", 'FontSize',14);
legend("\mu = "+ mu1, "\mu = " + mu2, 'FontSize',12);




