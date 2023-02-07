% Robust Regression
% OLS and PCR solutions for regression matrix estimation

clc;
clear;
close all;

rng(0);
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

%% Compare the estimation error
numOfTrials = 5000;
errorOLS_overall = zeros(numOfTrials,1);
errorPCR_overall = zeros(numOfTrials,1);
errorOLS = zeros(numOfTrials,size(Y,2));
errorPCR = zeros(numOfTrials,size(Y,2));

for i = 1:numOfTrials
    % get the test and estimated Y
    [Yols_test, Yols_estimate] = regval(Bols);
    [Ypcr_test, Ypcr_estimate] = regval(Bpcr);

    % compute the MSE error
    errorOLS_overall(i) = mean((Yols_test(:) - Yols_estimate(:)).^2);
    errorPCR_overall(i) = mean((Ypcr_test(:) - Ypcr_estimate(:)).^2);

    errorOLS(i,:) = mean((Yols_test - Yols_estimate).^2);
    errorPCR(i,:) = mean((Ypcr_test - Ypcr_estimate).^2);
end

%% display the comparison results
resultsOLS = [mean(errorOLS) mean(errorOLS_overall)];
resultsPCR = [mean(errorPCR) mean(errorPCR_overall)];

disp("The OLS estimation MSE error = " + mean(errorOLS_overall));
disp("The PCR estimation MSE error = " + mean(errorPCR_overall));
disp("The OLS error is " + ...
    100 * abs(mean(errorOLS_overall)-mean(errorPCR_overall)) / mean(errorOLS_overall) + ...
    "% larger than PCR error");

figure(1);
stem(1:length(mean(errorOLS)), mean(errorOLS), 'LineWidth',1); hold on;
stem(1:length(mean(errorPCR)), mean(errorPCR), 'X', 'LineWidth',1); hold off;
tit = title("MSE error between $\bf Y$ and $\bf \tilde{Y}_{test-OLS}$ and " + ...
    "$\bf \tilde{Y}_{test-PCR}$", 'FontSize', 20);
set(tit, 'Interpreter', 'latex');
leg = legend("MSE between $\bf Y$ and $\bf \tilde{Y}_{OLS}$", ...
    "MSE between $\bf Y$ and $\bf \tilde{Y}_{PCR}$", 'FontSize',16);
set(leg, 'Interpreter', 'latex');
xlabel("Column Index", 'FontSize',14); ylabel("Mean Square Error", 'FontSize',14);