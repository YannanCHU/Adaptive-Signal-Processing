clc;
clear;
close all;

load("EEG_Data_Assignment2.mat");

sPOz = POz;
t = (0:length(POz)-1)/fs;
N = length(POz);
figure(1); plot(t, POz); title("POz signal");
xlabel("Time (s)"); ylabel("EEG value");

% spectrogram parameters
windowSize = 5 * fs;
overlapPercent = 0.5;       % percentage of overlap
noverlap = round(overlapPercent * windowSize);
nFFT = 2^13;                % number of DFT points - determine the resolution of spectrogram

figure(2);
spectrogram(POz, windowSize, noverlap, nFFT, fs,'yaxis');
title("The spectrogram of a noisy single-channel EEG signal");
ylim([0,60]);

% generating a synthetic reference input composed of a sinusoid of 50 Hz corrupted
% by white Gaussian noise
noiseVar = 0.001;
noise = sqrt(noiseVar) * randn(1,N);
sine_clean = sin(2*pi*50*t);
sine_noisy = (sine_clean + noise).';
figure(3); plot(t, sine_noisy); title("a sinusoid of 50 Hz corrupted by white Gaussian noise");
xlabel("Time (s)"); ylabel("EEG value");

%% ANC parameters
filterOrders = [3,8,13,18];
stepSizes = [0.0001, 0.001, 0.01];
sPOZ_hat_ANC = cell(length(filterOrders), length(stepSizes));

for k = 1:1:length(filterOrders)
    M = filterOrders(k);
    for j = 1:1:length(stepSizes)
        mu = stepSizes(j);
        [noiseEst, ~, ~] = LMS_ale_anc(sine_noisy, sPOz, mu, M, [], 'ANC');
        sPOZ_hat_ANC{k,j} = sPOz - noiseEst;
    end
end

% display the results
figure(4);
tiledlayout(length(stepSizes), length(filterOrders),'TileSpacing','compact');
for j = 1:1:length(stepSizes)
    for k = 1:1:length(filterOrders)
        nexttile;
        spectrogram(sPOZ_hat_ANC{k,j}, windowSize, noverlap, nFFT, fs,'yaxis');
        title("The spectrogram of EEG signal (M = " +filterOrders(k)+", \mu = "+stepSizes(j)+")");
        ylim([0,60]);
    end
end

%% Apply the standard periodogram approach to the entire recording
% recommended frequency resolution -  5 DFT samples per Hz
f_res = 1/5;
[pxx_st,f_st] = periodogram(sPOz, hamming(length(sPOz)), fs/f_res, fs);
[pxx_st_anc,~] = periodogram(sPOZ_hat_ANC{2,2}, hamming(length(sPOz)), fs/f_res, fs);

figure(5);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(t, POz, t, sPOZ_hat_ANC{2,2}); title("Raw and Denoised POz signal", 'FontSize', 14);
xlabel("Time (s)", 'FontSize', 12); ylabel("EEG value", 'FontSize', 12);
legend("Raw", "Denoised", 'FontSize', 11);

nexttile;
plot(f_st, 10*log10(pxx_st) - 10*log10(pxx_st_anc), 'LineWidth', 2);
title('Difference between Periodograms of noisy and denoised EEG data', 'FontSize', 14);
xlabel("Frequency (Hz)", 'FontSize', 12);
ylabel('PSD difference (dB/Hz)', 'FontSize', 12);
xlim([0, 60]);
