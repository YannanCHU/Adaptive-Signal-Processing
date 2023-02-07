% Robust Regression
% singular values of X and Xnoise
% rank of X
% square error between singular values of X and Xnoise

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
r = 3;
SxnVec = diag(Sxn);
Sxn_lr = wthresh(Sxn, 'h', SxnVec(r+1));

% reconstruct the low rank approximated (denoised) matrix
Xdenoised = Uxn * Sxn_lr * Vxn';
% Xdenoised = Uxn(:,1:r) * Sxn(1:r,1:r) * Vxn(:,1:r)';
%% Mean Square Error Computation. 
% Error between the variables (columns) of the noiseless input matrix, X, 
% and those in the noise corrupted matrix Xnoise 
% and denoised matrix Xdenoised

errorX_Xnoise = mean((X-Xnoise).^2);
errorX_Xdenoised = mean((X-Xdenoised).^2);

figure();
stem(1:length(errorX_Xnoise), errorX_Xnoise, 'LineWidth', 4); hold on;
stem(1:length(errorX_Xdenoised), errorX_Xdenoised, 'LineWidth', 4); hold off;
tit = title("The Mean Square Error (MSE) between variables of noiseless $\bf X$ and those in" + ...
    "noisy $\bf X_{noise}$ and denoise $\bf \tilde{X}_{noise}$", 'FontSize',20);
set(tit, 'Interpreter', 'latex');
xlabel("Column Index", 'FontSize',16); ylabel("Mean Square Error", 'FontSize',16);
leg = legend("MSE between $\bf X$ and $\bf X_{noise}$", ...
    "MSE between $\bf X$ and $\bf \tilde{X}_{noise}$", 'FontSize',16);
set(leg, 'Interpreter', 'latex');