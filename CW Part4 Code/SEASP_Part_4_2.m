clc;
clear;
close all;

load("time-series.mat");

y_0_mean = y - mean(y);

N = length(y);      % the number of samples
mu = 1e-5;          % step size
AR_order = 4;       % AR(4)

[y_hat, errors, weights] = LMS_tanh([], y_0_mean, mu, AR_order, 0, 1, [], []);

startIndex = 1;
MSE = mean(errors(startIndex:end).^2);    % MSE
MSE_db = pow2db(MSE);

yhatVariance = var(y_hat(startIndex:end));      % variance of the output of LMS
yErrVariance = var(errors(startIndex:end));     % variance of the error
Rp = yhatVariance / yErrVariance;   % prediction gain
Rp_db = pow2db(Rp);

disp("The MSE is " + MSE_db);
disp("The Rp is " + Rp_db);



figure(1);
plot(1:N, y_0_mean, 1:N, y_hat, 'LineWidth', 2);
title(["Zero-mean version of the signal y and the one-step ahead prediction", ...
    sprintf("(MSE = %.3f dB, Rp = %.3f dB)", MSE_db, Rp_db) ], 'FontSize', 14);
xlabel("Sample Index", 'FontSize', 14);
ylabel("Signal amplitude", 'FontSize', 14);
legend("Zero-mean y", "One-step ahead prediction $\hat{y}$",'Interpreter','latex', 'FontSize', 14);
