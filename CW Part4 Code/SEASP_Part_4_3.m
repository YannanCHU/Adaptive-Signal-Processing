clc;
clear;
close all;

load("time-series.mat");

y_0_mean = y - mean(y);

N = length(y);      % the number of samples
mu = 1.5e-7;          % step size
AR_order = 4;       % AR(4)
startIndex = 1;

scalingFactors = 1:1:100;
MSEs_dB = zeros(length(scalingFactors),1);
Rps_db  = zeros(length(scalingFactors),1);

for i = 1:length(scalingFactors)
    [y_hat, errors, weights] = LMS_tanh([], y_0_mean, mu, AR_order, 0, scalingFactors(i), [], []);
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
legend("MSE", "Rp", 'FontSize', 14);

[minMSE, minMSEIdx] = min(MSEs_dB);
[maxRp, maxRpIdx] = max(Rps_db);

datatip(p1, scalingFactors(minMSEIdx), minMSE, 'Location', 'southwest');
datatip(p2, scalingFactors(maxRpIdx), maxRp, 'Location', 'northwest');

%%
[y_hat, errors, weights] = LMS_tanh([], y_0_mean, mu, AR_order, 0, scalingFactors(minMSEIdx), [], []);
MSE = mean(errors(startIndex:end).^2);    % MSE
MSE_db = pow2db(MSE);

yhatVariance = var(y_hat(startIndex:end));      % variance of the output of LMS
yErrVariance = var(errors(startIndex:end));     % variance of the error
Rp = yhatVariance / yErrVariance;   % prediction gain
Rp_db = pow2db(Rp);

disp("The overall MSE is " + MSE_db);
disp("The overall Rp is " + Rp_db);
disp("The optimal scaling factor - a is " + scalingFactors(minMSEIdx));

figure(2);
plot(1:N, y_0_mean, 1:N, y_hat, 'LineWidth', 2);
title(["Zero-mean version of the signal y and the one-step ahead prediction", ...
    sprintf("(MSE = %.3f dB, Rp = %.3f dB)", MSE_db, Rp_db) ], 'FontSize', 14);
xlabel("Sample Index", 'FontSize', 14);
ylabel("Signal amplitude", 'FontSize', 14);
legend("Zero-mean y", "One-step ahead prediction $\hat{y}$",'Interpreter','latex', 'FontSize', 14);

figure(3);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(1:N, pow2db(errors.^2), 'LineWidth', 2); hold on;
xline(180,'k--', 'LineWidth', 2); hold off;
xlabel("Sample Index", 'FontSize', 14);
ylabel("Squared Prediction Error in dB", 'FontSize', 14);
title("The squared prediction errors in each time step", 'FontSize', 14);

nexttile;
plot(1:N+1, weights, 'LineWidth', 2); hold on;
xline(180,'k--', 'LineWidth', 2); hold off;
xlabel("Sample Index", 'FontSize', 14);
ylabel("Weights", 'FontSize', 14);
title("The weight coefficients in each time step", 'FontSize', 14);
legend("w0", "w1", "w2", "w3", "w4", 'Orientation', 'horizontal', 'FontSize', 14);
xlim([1,N]);