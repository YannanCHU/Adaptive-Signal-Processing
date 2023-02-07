% adaptive filter coefficient estimathon based on leaky LMS

clc;
clear;
close all;

% number of samples
N = 1000;
% step sizes
mu1 = 0.05;
mu2 = 0.01;
mus = [mu1, mu2];
% leakage coefficient
gammas = 0.1:0.2:0.9;
% coefficients of AR model
a1 = 0.1;
a2 = 0.8;
% noise variance
noiseVar = 0.25;

% get the AR(2) model
% x(n) = 0 + a1*x(n-1) + a2*x(n-2) + noise
AR_MDL = arima('Constant',0,'AR',{a1 a2},'Variance',noiseVar);

%% 100 realisations and the averaged learning curve
% generate 100 realisations of x(n)
rng(0);
numOfRealisations = 100;
x100 = simulate(AR_MDL, N, 'NumPaths', numOfRealisations);   % simulate AR2 process
weightVecs = cell(length(mus), length(gammas), numOfRealisations);   % 2 x 5 x 100 cell storing the weight vector in last iteration

for j = 1:numOfRealisations
    for i = 1:length(mus)
        for k = 1:length(gammas)
            [~, errors, weights] = leaky_LMS_ar(x100(:,j), mus(i), gammas(k), 2);
            weightVecs{i,k,j} = weights(:, end);
        end
    end
end

referenceWeight = [a1;a2] .* ones(2,N);
averagedWeight = cell(length(mus), length(gammas));

for i = 1:length(mus)
    for k = 1:length(gammas)
        averagedWeight{i,k} = mean(cat(3,weightVecs{i,k,:}), 3);
    end
end

averagedWeight_mat = cell2mat(averagedWeight);
