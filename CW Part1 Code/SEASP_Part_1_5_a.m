% Real World Signals: Respiratory Sinus Arrhythmia from RR-Intervals
% Standard periodogram as well as the averaged periodogram with different window lengths
clc;
clear;
close all;

status = false;
if status
    load("./ECG_Data/RAW-ECG.mat");
    Trial_1 = data(770:2.5e5);
    Trial_2 = data(2.564e5:4.99e5);
    Trial_3 = data(5.228e5:7.473e5);
    [xRRI1,fsRRI1]=ECG_to_RRI(Trial_1, fs);
    print("T1 over");
    [xRRI2,fsRRI2]=ECG_to_RRI(Trial_2, fs);
    print("T2 over");
    [xRRI3,fsRRI3]=ECG_to_RRI(Trial_3, fs);
    print("T3 over");
end


%% import the RRI data
load("./ECG_Data/RRI-Data.mat");
x1 = detrend(xRRI1-mean(xRRI1));
x2 = detrend(xRRI2-mean(xRRI2));
x3 = detrend(xRRI3-mean(xRRI3));

% sampling frequency - prior knowledge
fs = 4;

figure(); plot((1:length(x1))/fs, x1, ...
    (1:length(x2))/fs, x2, (1:length(x3))/fs, x3, 'LineWidth',2);
title("Three different RRI signals");
xlabel("Time (s)"); ylabel("Data");
legend("RRI1", "RRI2", "RRI3");

%% standard periodogram
nfft = 2048;
[Pxx1_st,f_st] = periodogram(x1, hamming(length(x1)), nfft, fs);
[Pxx2_st,~] = periodogram(x2, hamming(length(x2)), nfft, fs);
[Pxx3_st,~] = periodogram(x3, hamming(length(x3)), nfft, fs);

%% average periodogram
windowSizes = [50, 150];
[Pxx1_avg_w50, ~] = pwelch(x1, windowSizes(1)*fs, 0, nfft, fs);
[Pxx1_avg_w150, ~] = pwelch(x1, windowSizes(2)*fs, 0, nfft, fs);

[Pxx2_avg_w50, ~] = pwelch(x2, windowSizes(1)*fs, 0, nfft, fs);
[Pxx2_avg_w150, ~] = pwelch(x2, windowSizes(2)*fs, 0, nfft, fs);

[Pxx3_avg_w50, ~] = pwelch(x3, windowSizes(1)*fs, 0, nfft, fs);
[Pxx3_avg_w150, ~] = pwelch(x3, windowSizes(2)*fs, 0, nfft, fs);

%% display the PSD estimates
figure();
subplot(311);
plot(f_st, pow2db(Pxx1_st), f_st, pow2db(Pxx1_avg_w50), ...
    f_st, pow2db(Pxx1_avg_w150), 'LineWidth', 2);
title("Standard and Averaged Periodogram of RRI data (Trial 1)",'FontSize',13);
xlabel("Frequency (Hz)",'FontSize',13);
ylabel("Power/Frequency (dB/Hz)",'FontSize',13);
legend("Standard", "Averaged (window = 50)", "Averaged (window = 150)",'FontSize',11);

subplot(312);
plot(f_st, pow2db(Pxx2_st), f_st, pow2db(Pxx2_avg_w50), ...
    f_st, pow2db(Pxx2_avg_w150), 'LineWidth', 2);
title("Standard and Averaged Periodogram of RRI data (Trial 2)",'FontSize',13);
xlabel("Frequency (Hz)",'FontSize',13);
ylabel("Power/Frequency (dB/Hz)",'FontSize',13);
legend("Standard", "Averaged (window = 50)", "Averaged (window = 150)",'FontSize',11);

subplot(313);
plot(f_st, pow2db(Pxx3_st), f_st, pow2db(Pxx3_avg_w50), ...
    f_st, pow2db(Pxx3_avg_w150), 'LineWidth', 2);
title("Standard and Averaged Periodogram of RRI data (Trial 3)",'FontSize',13);
xlabel("Frequency (Hz)",'FontSize',13);
ylabel("Power/Frequency (dB/Hz)",'FontSize',13);
legend("Standard", "Averaged (window = 50)", "Averaged (window = 150)",'FontSize',11);

