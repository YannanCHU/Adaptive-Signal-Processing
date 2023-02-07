% Real World Signals: Respiratory Sinus Arrhythmia from RR-Intervals
% AR spectrum estimate
clc;
clear;
close all;

%% import the RRI data
load("./ECG_Data/RRI-Data.mat");
x1 = detrend(xRRI1-mean(xRRI1));
x2 = detrend(xRRI2-mean(xRRI2));
x3 = detrend(xRRI3-mean(xRRI3));

% sampling frequency - prior knowledge
fs = 4;
% number of frequency bins
nfft = 2048;

%% standard periodogram
[Pxx1_st,f_st] = periodogram(x1, hamming(length(x1)), nfft, fs);
[Pxx2_st,~] = periodogram(x2, hamming(length(x2)), nfft, fs);
[Pxx3_st,~] = periodogram(x3, hamming(length(x3)), nfft, fs);

%% AR spectrum estimation
% order of underlying AR model
p = 2:1:36;
% create three cell to store estimated PSD
PxxAREst1 = cell(length(p),1);
PxxAREst2 = cell(length(p),1);
PxxAREst3 = cell(length(p),1);

for i = 1:1:length(p)
    % normalized autoregressive (AR) parameters corresponding to a model of
    % order p for the input array x;
    [ARparameters1,noiseVarEst1,~] = aryule(x1,p(i));
    % estimated PSD of AR model
    [he1,f] = freqz(sqrt(noiseVarEst1),ARparameters1,nfft,fs);   % get the frequency reeponse of filter
    PxxAREst1{i} = pow2db(abs(he1).^2);        % get the true PSD

    % normalized autoregressive (AR) parameters corresponding to a model of
    % order p for the input array x;
    [ARparameters2,noiseVarEst2,~] = aryule(x2,p(i));
    % estimated PSD of AR model
    he2 = freqz(sqrt(noiseVarEst2),ARparameters2,nfft,fs);   % get the frequency reeponse of filter
    PxxAREst2{i} = pow2db(abs(he2).^2);        % get the true PSD

    % normalized autoregressive (AR) parameters corresponding to a model of
    % order p for the input array x;
    [ARparameters3,noiseVarEst3,~] = aryule(x3,p(i));
    % estimated PSD of AR model
    he3 = freqz(sqrt(noiseVarEst3),ARparameters3,nfft,fs);   % get the frequency reeponse of filter
    PxxAREst3{i} = pow2db(abs(he3).^2);        % get the true PSD
end

%% Display the estimation results
index1 = 23;
index2 = 12;
index3 = 7;

figure(1);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
plot(f_st, pow2db(Pxx1_st), f,PxxAREst1{index1}, f,PxxAREst1{index2}, f,PxxAREst1{index3}, 'LineWidth', 2);
title("AR spectrum estimate of RRI data (trial 1)");
xlabel("Frequency (Hz)"); ylabel("Power/Frequency (dB/Hz)");
legend("Periodogram", "Order = " + p(index1), "Order = " + p(index2), "Order = " + p(index3));

nexttile;
plot(f_st, pow2db(Pxx2_st), f,PxxAREst2{index1}, f,PxxAREst2{index2}, f,PxxAREst2{index3}, 'LineWidth', 2);
title("AR spectrum estimate of RRI data (trial 2)");
xlabel("Frequency (Hz)"); ylabel("Power/Frequency (dB/Hz)");
legend("Periodogram", "Order = " + p(index1), "Order = " + p(index2), "Order = " + p(index3));

nexttile;
plot(f_st, pow2db(Pxx3_st), f,PxxAREst3{index1}, f,PxxAREst3{index2}, f,PxxAREst3{index3}, 'LineWidth', 2);
title("AR spectrum estimate of RRI data (trial 3)");
xlabel("Frequency (Hz)"); ylabel("Power/Frequency (dB/Hz)");
legend("Periodogram", "Order = " + p(index1), "Order = " + p(index2), "Order = " + p(index3));

figure(2);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
plot(f_st, pow2db(Pxx1_st), f,PxxAREst1{index1}, 'LineWidth', 2);
title("AR spectrum estimate of RRI data (trial 1)",'FontSize',13);
xlabel("Frequency (Hz)",'FontSize',13); ylabel("Power/Frequency (dB/Hz)",'FontSize',13);
legend("Periodogram", "AR (" + p(index1) + ") Spectrum",'FontSize',11);

nexttile;
plot(f_st, pow2db(Pxx2_st), f,PxxAREst2{index2}, 'LineWidth', 2);
title("AR spectrum estimate of RRI data (trial 2)",'FontSize',13);
xlabel("Frequency (Hz)",'FontSize',13); ylabel("Power/Frequency (dB/Hz)",'FontSize',13);
legend("Periodogram", "AR (" + p(index2) + ") Spectrum",'FontSize',11);

nexttile;
plot(f_st, pow2db(Pxx3_st), f,PxxAREst3{index3}, 'LineWidth', 2);
title("AR spectrum estimate of RRI data (trial 3)",'FontSize',13);
xlabel("Frequency (Hz)",'FontSize',13); ylabel("Power/Frequency (dB/Hz)",'FontSize',13);
legend("Periodogram", "AR (" + p(index3) + ") Spectrum",'FontSize',11);
