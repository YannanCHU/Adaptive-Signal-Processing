% 1.2 Periodogram-based Methods Applied to Real–World Data
% b) an electroencephalogram (EEG) experiment
clc;
clear;
close all;

%% import the EEG data
load("./EEG_Data/EEG_Data_Assignment1.mat");
x = POz;
figure(1); plot([0:length(POz)-1]/fs, POz); title("POz signal");
xlabel("Time (s)"); ylabel("EEG value");
% number of samples
numOfSample = length(x);
% recommended frequency resolution -  5 DFT samples per Hz
f_res = 1/5;

%% Apply the standard periodogram approach to the entire recording
[pxx_st,f_st] = periodogram(x, hamming(length(x)), fs/f_res, fs);

figure(2);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
p1 = plot(f_st, 10*log10(pxx_st), 'LineWidth', 2);
xlabel("Frequency (Hz)",'FontSize',12);
ylabel('Power/frequency (dB/(rad/sample))','FontSize',12);
title('Standard Periodogram of EEG data','FontSize',13);
xlim([0, 60]);
for i = 1:5
    datatip(p1,'DataIndex',1+13*5*i,'Location','northeast');
end
%% Averaged periodogram
% specify three different window sizes (hamming window is used)
windowSizes = [1, 5, 10]; % unit: second
% create an 2D array to store the results
pxx_avg_array = zeros(0.5*fs/f_res+1, length(windowSizes));

for i = 1:length(windowSizes)
    [pxx_avg_array(:,i), ~] = pwelch(x, windowSizes(i)*fs, 0, fs/f_res, fs);
end

% figure(2);
nexttile;
p2 = plot(f_st, 10*log10(pxx_avg_array), 'LineWidth', 2);
legend("Window length = 1 s", "Window length = 5 s", "Window length = 10 s",'FontSize',11)
xlabel("Frequency (Hz)",'FontSize',12);
ylabel('PSD (dB/(rad/sample))','FontSize',12);
title('Averaged periodogram of EEG data with different window lengths','FontSize',13);
xlim([0, 60]);
for i = 1:5
    datatip(p2(3,1), "DataIndex", 1+13*5*i, 'Location', 'northeast');
%     datatip(p2(1,1), 'DataIndex', 1+13*5*i, 'Location', 'southeast');
end