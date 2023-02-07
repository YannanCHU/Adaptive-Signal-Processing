clc;
clear;
close all;

%% generate the frequency modulated (FM) signal - yn
noiseVar = 0.05;    % noise variance
N = 1500;           % number of signals
n = 1:1:N;          % the sequence of smaple indexes
fs = 1500;          % the sampling frequency
fn = (1<=n) .* (n<=500) .* 100 + ...
     (501<=n) .* (n<=1000) .* (100 + (n-500)/2) + ...
     (1001<=n) .* (n<=1500) .* (100 + ((n-1000)/25).^2);
Phi = cumsum(fn);   % the phase signal

rng(12);
eta = sqrt(noiseVar / 2) * (randn(1,N) + 1j * randn(1,N));    % the noise signal
coefficient = circularityMeasure(eta);
disp("The circularity of complex noise is " + coefficient);
yn = exp(1j*( 2*pi*Phi(n) / fs )) + eta;    % the FM signal

%% AR based method
AR_order = 1;       % the order of AR filter
mus = [0.03];      % step size
numOfPoints = 1024;
CLMS_sepctrums = cell(length(mus),1);

for k = 1:1:length(mus)
    H = zeros(numOfPoints, N);
    [~, ~, AR_coefs] = CLMS_arma([], yn, mus(k), AR_order, 0);
    for i = 1:N
        %  Run complex-valued LMS algorithm to estimate AR coefficient
        [hFreResp, wAngFreResp] = freqz(1, [1; -conj(AR_coefs(i))], numOfPoints, fs);  % compute power spectrum
        H(:,i) = abs(hFreResp) .^ 2;    % store it in a matrix
    end

    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    CLMS_sepctrums{k,1} = H;
end

%% DFT based method
mu_dft = 1;
gammas = [0, 0.1];
DFT_spectrum = cell(length(gammas), 1);

for i = 1:1:length(gammas)
    DFT_H = zeros(N, N);
    [~, ~, DFT_coefs] = DFT_CLMS(yn, mu_dft, gammas(i));
    DFT_H = abs(DFT_coefs).^2;
    % Remove outliers in the matrix H
    medianDFT_H = 50*median(median(DFT_H));
    DFT_H(DFT_H > medianDFT_H) = medianDFT_H;
    DFT_spectrum{i,1} = DFT_H;
end

figure();
tiledlayout(2,1,'TileSpacing','compact');

for i = 1:1:length(gammas)
    nexttile;
    surf(n, linspace(0,fs,N), DFT_spectrum{i,1},'LineStyle','none');     view(2);
    title("DFT-CLMS time-frequency spectrum estimate (\mu = " + mu_dft +", \gamma = " + gammas(i) + ")", 'FontSize',14);
    c = colorbar;   ylabel(c, "PSD", 'FontSize',12);
    xlabel("Sample Index", 'FontSize',14);
    ylabel("Frequency (Hz)", 'FontSize',14);
    ylim([0, fs/2]);
end