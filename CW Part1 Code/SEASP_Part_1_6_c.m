% Robust Regression
% OLS and PCR solutions for regression matrix estimation

clc;
clear;
close all;

%% import the PCR data
load("./PCR/PCAPCR");

%% Singular Value Decomposition
% SVD of X
[Ux,Sx,Vx] = svd(X, 'econ');
% SVD of Xnoise
[Uxn,Sxn,Vxn] = svd(Xnoise, 'econ');

%% create a low-rank approximation of Xnoise
% remove the singular values corresponding to noise subspaces
% reconstruct the low rank approximated (denoised) matrix
r = 3;
Xdenoised = Uxn(:,1:r) * Sxn(1:r,1:r) * Vxn(:,1:r)';

%% OLS and PCR solutions - training set 
Bols = (Xnoise.' * Xnoise) \ Xnoise.' * Y;
Bpcr = Vxn(:,1:r) * inv(Sxn(1:r,1:r)) * Uxn(:,1:r).' * Y;
% Bpcr = Vxn(:,1:r) * (Sxn(1:r,1:r) \ Uxn(:,1:r).') * Y;

%% Compare the estimation error - training set 
Yols_training = Xnoise * Bols;
Ypcr_training = Xdenoised * Bpcr;

errorOLS_training = mean((Y-Yols_training).^2);
errorOLS_training_overall = mean((Y(:)-Yols_training(:)).^2);
errorPCR_training = mean((Y-Ypcr_training).^2);
errorPCR_training_overall = mean((Y(:)-Ypcr_training(:)).^2);

figure();
stem(1:length(errorOLS_training), errorOLS_training, 'LineWidth',1); hold on;
stem(1:length(errorPCR_training), errorPCR_training, 'X', 'LineWidth',1); hold off;
tit = title("MSE error between $\bf Y$ and $\bf \tilde{Y}_{OLS}$ and " + ...
    "$\bf \tilde{Y}_{PCR}$", 'FontSize', 20);
set(tit, 'Interpreter', 'latex');
leg = legend("MSE between $\bf Y$ and $\bf \tilde{Y}_{OLS}$", ...
    "MSE between $\bf Y$ and $\bf \tilde{Y}_{PCR}$", 'FontSize',16);
set(leg, 'Interpreter', 'latex');
xlabel("Column Index", 'FontSize',14); ylabel("Mean Square Error", 'FontSize',14);

%% Compare the estimation error - test set
% estimated output
Yols_test = Xtest * Bols;
Ypcr_test = Xtest * Bpcr;

errorOLS_test = mean((Ytest-Yols_test).^2);
errorOLS_test_overall = mean((Ytest(:)-Yols_test(:)).^2);
errorPCR_test = mean((Ytest-Ypcr_test).^2);
errorPCR_test_overall = mean((Ytest(:)-Ypcr_test(:)).^2);

disp("PCR error is " + 100*abs(errorPCR_test_overall-errorOLS_test_overall)/errorOLS_test_overall + ...
    "% better than OLS error");

figure();
stem(1:length(errorOLS_test), errorOLS_test, 'LineWidth',1); hold on;
stem(1:length(errorPCR_test), errorPCR_test, 'X', 'LineWidth',1); hold off;
tit = title("MSE error between $\bf Y$ and $\bf \tilde{Y}_{test-OLS}$ and " + ...
    "$\bf \tilde{Y}_{test-PCR}$", 'FontSize', 20);
set(tit, 'Interpreter', 'latex');
leg = legend("MSE between $\bf Y$ and $\bf \tilde{Y}_{OLS}$", ...
    "MSE between $\bf Y$ and $\bf \tilde{Y}_{PCR}$", 'FontSize',16);
set(leg, 'Interpreter', 'latex');
xlabel("Column Index", 'FontSize',14); ylabel("Mean Square Error", 'FontSize',14);
