clc;
clear;
close all;

load("time-series.mat");

N = length(y);      % the number of samples
mu = 1.5e-7;          % step size
AR_order = 4;       % AR(4)
startIndex = 1;

scalingFactors = 1:1:100;
MSEs_dB = zeros(length(scalingFactors),1);
Rps_db  = zeros(length(scalingFactors),1);

for i = 1:length(scalingFactors)
    [y_hat, errors, ~] = LMS_tanh([], y, mu, AR_order, 0, scalingFactors(i), 'bias', []);
    MSEs_dB(i) = pow2db(mean(errors(startIndex:end).^2));    % MSE
    
    yhatVariance = var(y_hat(startIndex:end));      % variance of the output of LMS
    yErrVariance = var(errors(startIndex:end));     % variance of the error
    Rps_db(i) = pow2db(yhatVariance / yErrVariance);      % prediction gain
end

figure();
p1 = plot(scalingFactors, MSEs_dB, 'LineWidth', 2); hold on;
p2 = plot(scalingFactors, Rps_db,  'LineWidth', 2); hold off;
xlabel("Sacling Factor - a", 'FontSize', 14);
ylabel("MSE (blue) in dB and Rp (red) in dB", 'FontSize', 14);
title("The Mean Square Error (MSE) and Prediction gain (Rp) against Sacling Factor (a)", 'FontSize', 14);
legend("MSE", "Rp");

[minMSE, minMSEIdx] = min(MSEs_dB);
[maxRp, maxRpIdx] = max(Rps_db);

datatip(p1, scalingFactors(minMSEIdx), minMSE, 'Location', 'southwest');
datatip(p2, scalingFactors(maxRpIdx), maxRp);

a_opt = scalingFactors(minMSEIdx);

%% pretrain the weights by using 100 iterations and 20 samples
weights_init = zeros(5,1);
numOfIters = 100;
numOfSamples = 20;
for k = 1:numOfIters
    [~, ~, weights] = LMS_tanh([], y(1:numOfSamples), mu, AR_order, 0, scalingFactors(minMSEIdx), 'bias', weights_init);
    weights_init = weights(:,end);
end

% the MSE error produced by following method is lower than above method
% pretrainingSamples = repmat(y(1:numOfSamples), numOfIters, 1);
% [~, ~, weights2] = LMS_tanh([], pretrainingSamples, mu, AR_order, 0, scalingFactors(minMSEIdx), 'bias', zeros(5,1));
% weights_init = weights2(:,end);

%% predict the complete signal
[y_hat, errors, weights] = LMS_tanh([], y, mu, AR_order, 0, scalingFactors(minMSEIdx), 'bias', weights_init);
MSE = mean(errors(startIndex:end).^2);    % MSE
MSE_db = pow2db(MSE);

yhatVariance = var(y_hat(startIndex:end));      % variance of the output of LMS
yErrVariance = var(errors(startIndex:end));     % variance of the error
Rp = yhatVariance / yErrVariance;   % prediction gain
Rp_db = pow2db(Rp);

disp("The MSE is " + MSE_db);
disp("The Rp is " + Rp_db);

figure(2);
plot(1:N, y, 1:N, y_hat, 'LineWidth', 2);
title([sprintf("Signal y and the one-step ahead prediction (Biased and scaled (a = %d) tanh LMS with pre-trained initial weights)", a_opt), ...
    sprintf("(MSE = %.3f dB, Rp = %.3f dB)", MSE_db, Rp_db) ], 'FontSize', 14);
xlabel("Sample Index", 'FontSize', 14);
ylabel("Signal amplitude", 'FontSize', 14);
legend("Zero-mean y", "One-step ahead prediction $\hat{y}$",'Interpreter','latex', 'FontSize', 14);
