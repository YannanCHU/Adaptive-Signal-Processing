clc;
clear;
close all;

load("time-series.mat");

y_0_mean = y - mean(y);

N = length(y);      % the number of samples
mu = 1e-5;          % step size
AR_order = 4;       % AR(4)

[y_hat, errors, weights] = LMS_arma([], y_0_mean, mu, AR_order, 0);

steadyStartIndex = 1;
MSE = mean(errors(steadyStartIndex:end).^2);    % MSE
MSE_db = pow2db(MSE);

yhatVariance = var(y_hat(steadyStartIndex:end));      % variance of the output of LMS
yErrVariance = var(errors(steadyStartIndex:end));     % variance of the error
Rp = yhatVariance / yErrVariance;   % prediction gain
Rp_db = pow2db(Rp);

disp("The overall MSE is " + MSE_db);
disp("The overall Rp is " + Rp_db);
disp("The initial MSE (1-200) is " + pow2db(mean(errors(1:200).^2)));
disp("The initial Rp (1-200) is " + pow2db(var(y_hat(1:200)) / var(errors(1:200))));
disp("The initial MSE (201-1000) is " + pow2db(mean(errors(201:1000).^2)));
disp("The initial Rp (201-1000) is " + pow2db(var(y_hat(201:1000)) / var(errors(201:1000))));

figure(1);
plot(1:N, y_0_mean, 1:N, y_hat, 'LineWidth', 2);
title(["Zero-mean version of the signal y and the one-step ahead prediction", ...
    sprintf("(MSE = %.3f dB, Rp = %.3f dB)", MSE_db, Rp_db) ], 'FontSize', 14);
xlabel("Sample Index", 'FontSize', 14);
ylabel("Signal amplitude", 'FontSize', 14);
legend("Zero-mean y", "One-step ahead prediction $\hat{y}$",'Interpreter','latex', 'FontSize', 14);
