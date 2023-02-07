% generate the PSD estimate of several realisations of a random process and plot them
% Generate different signals composed of sinusoids corrupted by noise and 
% elaborate on how disperse are the different realisation of the spectral estimate

clc;
clear;
close all;

seed = 10;
rng(seed);

%% Signal simulation
Fs = 80;                % sampling frequency
Ts = 1 / Fs;            % sampling interval
Ns = 2048;              % number of samples
ts = Ts * (0:1:Ns-1);   % sampling time

% Noiseless sinusoidal signals
fx1 = 12;       fx2 = 24;                % signal frequency
amp1 = 0.4;     amp2 = 0.7;              % signal amplitude
ph1 = 0;        ph2 = 0;                 % signal phase
noiselessSine = amp1 * sin(2*pi*fx1*ts+ph1) + amp2 * sin(2*pi*fx2*ts+ph2);

% The number of realisations
numOfRealisations = 500;
% The noise
noise = randn(Ns, numOfRealisations);
% The noisy sinusoidal signals
noisySines = noiselessSine + noise.';

% Creat matrices to store calculated ACF and PSD of each realisation
Rxxs = zeros(numOfRealisations, 2*Ns-1);
Pxxs = zeros(numOfRealisations, 2*Ns-1);

% Generate independent realisations of noisy sinusoidal signals
for i = 1:1:numOfRealisations
    % biased ACF estimator
    Rxxs(i,:) = xcorr(noisySines(i,:), "biased");
    % estimate PSD
    Pxxs(i,:) = abs(fftshift(fft(ifftshift(Rxxs(i,:)))));
end

% mean and standard deviation of the estimated PSDs
meanPxx = mean(Pxxs);
stdPxx = std(Pxxs);

%% display the results
% get the normalised freuqency
fs = Fs * (-Ns+1:1:Ns-1) / (2 * Ns);

figure(1);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;

p1(1:numOfRealisations) = plot(fs, Pxxs,"c", 'LineWidth', 2); hold on;
p1(numOfRealisations+1) = plot(fs, meanPxx, 'b', 'LineWidth', 2); hold off;
xlim([fs(Ns), fs(end)]);
title("PSD estimates (different realisations and mean)",'FontSize',14);
xlabel("Frequency (Hz)",'FontSize',14);
ylabel("PSD",'FontSize',14);
legend(p1([1,numOfRealisations+1]), 'realisations', 'mean','FontSize',11);

nexttile;
plot(fs, stdPxx,"r", 'LineWidth', 2);
xlim([fs(Ns), fs(end)]);
title("Standard deviation of the PSD estimate",'FontSize',14);
xlabel("Frequency (Hz)",'FontSize',14);
ylabel("PSD",'FontSize',14);

%% display the PSD in the unit of dB
figure(2);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
p2(1:numOfRealisations) = plot(fs, 10*log10(Pxxs),"c"); hold on;
p2(numOfRealisations+1) = plot(fs, mean(10*log10(Pxxs)), 'b', 'LineWidth', 2); hold off;
xlim([fs(Ns), fs(end)]);
title("PSD estimates (different realisations and mean)",'FontSize',14);
xlabel("Frequency (Hz)",'FontSize',14);
ylabel("PSD (dB)",'FontSize',14);
legend(p2([1,numOfRealisations+1]), 'realisations', 'mean','FontSize',11);

nexttile;
plot(fs, std(10*log10(Pxxs)),"r", 'LineWidth', 2);
xlim([fs(Ns), fs(end)]);
title("Standard deviation of the PSD estimate",'FontSize',14);
xlabel("Frequency (Hz)",'FontSize',14);
ylabel("PSD (dB)",'FontSize',14);
