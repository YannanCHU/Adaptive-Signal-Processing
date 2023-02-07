% steady state values of the adaptive filter coefficients

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
% set a time index corresponding to the steady state
% this value is delibrately set to be large
steadyTimeIdex = 500;

% get the AR(2) model
% x(n) = 0 + a1*x(n-1) + a2*x(n-2) + noise
AR_MDL = arima('Constant',0,'AR',{a1 a2},'Variance',noiseVar);

%% 100 realisations and the averaged learning curve
% generate 100 realisations of x(n)
rng(0);
numOfRealisations = 100;
x100 = simulate(AR_MDL, N, 'NumPaths', numOfRealisations);   % simulate AR2 process
weightVecs = cell(length(mus), numOfRealisations);   % 2 x 100 cell storing the weight vector in each iteration

for j = 1:numOfRealisations
    for i = 1:length(mus)
        [~, ~, weightVecs{i,j}] = LMS_arma([], x100(:,j), mus(i), 2, 0);
    end
end

referenceWeight = [a1;a2] .* ones(2,N);
averagedWeight_mu1 = mean(cat(3,weightVecs{1,:}),3);
averagedWeight_mu2 = mean(cat(3,weightVecs{2,:}),3);

figure();
plot(1:N, averagedWeight_mu1(1,:), 'r--', ...
    1:N, averagedWeight_mu2(1,:), 'g--', ...
    1:N, referenceWeight(1,:), 'b--', ...
    1:N, averagedWeight_mu1(2,:), 'r-', ...
    1:N, averagedWeight_mu2(2,:), 'g-', ...
    1:N, referenceWeight(2,:), 'b-', ...
    'LineWidth', 2);
ylim([0,0.9]);
title("Averaged weight coefficients over 100 trials", 'FontSize',16);
xlabel("Time index", 'FontSize',14);
ylabel("Weight coefficient", 'FontSize', 14);
legend("$\hat{a}_1 (\mu = 0.05)$", "$\hat{a}_1 (\mu = 0.01)$", "$a_1 (true)$", ...
    "$\hat{a}_2 (\mu = 0.05)$", "$\hat{a}_2 (\mu = 0.01)$", "$a_2 (true)$", ...
    'Interpreter','latex', 'Location', 'east','FontSize',14);

error_mu1 = 100*abs(averagedWeight_mu1(:,end)-[a1;a2]) ./ [a1;a2];
error_mu2 = 100*abs(averagedWeight_mu2(:,end)-[a1;a2]) ./ [a1;a2];

fprintf("When mu is 0.05, the coefficients are %.3f and %.3f\n", averagedWeight_mu1(:,end));
fprintf("The error rates are %.3f %% and %.3f %%\n\n", error_mu1);

fprintf("When mu is 0.01, the coefficients are %.3f and %.3f\n", averagedWeight_mu2(:,end));
fprintf("The error rates are %.3f %% and %.3f %%\n\n", error_mu2);