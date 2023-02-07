% misadjustment by time-averaging over the steady state of the
% ensemble-averaged learning curves using 100 independent trials

clc;
clear;
% close all;

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
MSEs = zeros(length(mus), numOfRealisations);   % 2 x 100 matrix storing the MSE of each trial

for j = 1:numOfRealisations
    for i = 1:length(mus)
        [~, errors, ~] = LMS_arma([], x100(:,j), mus(i), 2, 0);
        MSEs(i,j) = mean(errors(steadyTimeIdex:end) .^ 2);
    end
end

EMSE = mean(MSEs - noiseVar, 2);
misadjustments = EMSE ./ noiseVar;

%% misadjustment approximation of the LMS
A = [1, -a1, -a2; a1, a2-1, 0; a2, a1, -1];
b = [noiseVar; 0; 0];
rs = A \ b;

R = [rs(1), rs(2); rs(2), rs(1)];
misadjustments_lms = sum(diag(R)) * mus / 2;

misadjustments_diff = abs(misadjustments.' - misadjustments_lms);
misadjustments_diff_percent = 100*misadjustments_diff ./ misadjustments_lms;
disp(misadjustments_diff_percent + "%");

fprintf("When mu = %.2f, the MSE = %.4f, the EMSE = %.4f, the misadjustement = %.4f, the approximated misadjustement = %.4f;\n", ...
    mus(1), mean(MSEs(1,:)), EMSE(1), misadjustments(1), misadjustments_lms(1));
fprintf("The difference between estimated and theoretical misadjustement = %.4f \n\n", ...
    misadjustments_diff(1));
fprintf("When mu = %.2f, the MSE = %.4f, the EMSE = %.4f, the misadjustement = %.4f, the approximated misadjustement = %.4f;\n", ...
    mus(2), mean(MSEs(2,:)), EMSE(2), misadjustments(2), misadjustments_lms(2));
fprintf("The difference between estimated and theoretical misadjustement = %.4f \n\n", ...
    misadjustments_diff(2));

