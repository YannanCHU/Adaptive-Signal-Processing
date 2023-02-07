clc;
clear;
close all;

%% import the EEG data
load("EEG_Data_Assignment1.mat");
figure(1); plot([0:length(POz)-1]/fs, POz); title("POz signal");
xlabel("Time (s)"); ylabel("EEG value");

a = 2000;
yn = POz(a:a+1200-1);
N = 1200;

%% DFT based method
mu_dft = 1;
gammas = [0, 0.01];
DFT_spectrum = cell(length(gammas), 1);

for i = 1:1:length(gammas)
    DFT_H = zeros(N, N);
    [~, ~, DFT_coefs] = DFT_CLMS(yn.', mu_dft, gammas(i));
    DFT_H = abs(DFT_coefs).^2;
    % Remove outliers in the matrix H
    medianDFT_H = 300*median(median(DFT_H));
    DFT_H(DFT_H > medianDFT_H) = medianDFT_H;
    DFT_spectrum{i,1} = DFT_H;
end

figure();
tiledlayout(1,2,'TileSpacing','compact');

for i = 1:1:length(gammas)
    nexttile;
    surf(a:a+1200-1, linspace(0,fs,N), DFT_spectrum{i,1},'LineStyle','none');     view(2);
    title("DFT-CLMS time-frequency spectrum estimate (\mu = " + mu_dft +", \gamma = " + gammas(i) + ")", 'FontSize',14);    
    c = colorbar;   ylabel(c, "PSD", 'FontSize',12);
    xlabel("Sample Index", 'FontSize',14);
    ylabel("Frequency (Hz)", 'FontSize',14);
    ylim([0, 100]);
end
