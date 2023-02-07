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

%% Square Error between Sx and Sxn
errorSq = (diag(Sx) - diag(Sxn)).^2;

figure(1);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
stem(1:length(diag(Sx)), diag(Sx), 'LineWidth', 2); hold on;
stem(1:length(diag(Sxn)), diag(Sxn), 'x', 'LineWidth', 2); hold off;
xlabel("Index",'FontSize',13); ylabel("Singular Value",'FontSize',13);
title("Singular values of X and Xnoise",'FontSize',14);
legend("X", "Xnoise",'FontSize',12);

nexttile;
stem(1:length(errorSq), errorSq, 'LineWidth', 2); hold on;
xlabel("Index",'FontSize',13); ylabel("Square Error",'FontSize',13);
title("Square error between each singular value of X and Xnoise",'FontSize',14);

