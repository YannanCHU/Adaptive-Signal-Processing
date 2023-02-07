% 1.2 Periodogram-based Methods Applied to Real–World Data
% a) The sunspot time series
clc;
clear;
close all;

%% import the sunspot time series data
load sunspot.dat
ts = sunspot(:,1);              % sampling time
x_original = sunspot(:,2);      % original data
Ns = length(x_original);        % number of samples

% remove the mean and detrend the original data
x_unbiased_detrend = detrend(x_original - mean(x_original), 1);

% remove the mean from the logarithmic data 
epsilon = eps;    % a small value used to avoid log(0);
x_log_unbiased = log(x_original+epsilon) - mean(log(x_original+epsilon));

% plot the original and adjusted time series
figure(1);
subplot(311);
plot(ts, x_original, ts, x_unbiased_detrend, ts, x_log_unbiased, 'LineWidth', 2);
title("Sunspot time series",'FontSize',13);
legend("Original data", "remove mean and trend", "remove mean from log data",'FontSize',11);
xlabel("Cycle (year)");
ylabel("Number of sunspots");
xlim([min(ts), max(ts)]);

%% apply periodogram analysis to data without window
[pxx_original,f_original] = periodogram(x_original,[],[],1);
[pxx_unbiased_detrend,~] = periodogram(x_unbiased_detrend,[],[],1);
[pxx_log_unbiased,~] = periodogram(x_log_unbiased,[],[],1);

% plot the periodogram
subplot(312);
plot(f_original, 10*log10(pxx_original), ...
     f_original, 10*log10(pxx_unbiased_detrend), '--', ...
     f_original, 10*log10(pxx_log_unbiased), 'LineWidth', 2)
xlabel("Normalized Frequency  (\times\pi rad/sample)",'FontSize',12);
ylabel('PSD (dB/(rad/sample))','FontSize',12);
title('Periodogram (without window) of Relative Sunspot Number Data','FontSize',13);
legend("Original data", "remove mean and trend", "remove mean from log data",'FontSize',11);
ylim([-20, 60]);

%% specify a hamming window when applying periodogram analysis to data
[pxx_w_original,f_w_original] = periodogram(x_original,hamming(Ns),Ns,1);
[pxx_w_unbiased_detrend,~] = periodogram(x_unbiased_detrend,hamming(Ns),Ns,1);
[pxx_w_log_unbiased,~] = periodogram(x_log_unbiased,hamming(Ns),Ns,1);

% plot the periodogram
subplot(313);
plot(f_w_original, 10*log10(pxx_w_original), ...
     f_w_original, 10*log10(pxx_w_unbiased_detrend), '--', ...
     f_w_original, 10*log10(pxx_w_log_unbiased), 'LineWidth', 2)
xlabel("Normalized Frequency  (\times\pi rad/sample)",'FontSize',12);
ylabel('PSD (dB/(rad/sample))','FontSize',12);
title('Periodogram (with Hamming window) of Relative Sunspot Number Data','FontSize',13);
legend("Original data", "remove mean and trend", "remove mean from log data",'FontSize',11);

figure();
subplot(311);
plot(f_w_original, pxx_w_original);
xlabel("Normalized Frequency  (\times\pi rad/sample)");
ylabel('Power/frequency (W/(rad/sample))');
title("Periodogram (with Hamming window) of Original signal");
subplot(312);
plot(f_w_original, pxx_w_unbiased_detrend);
xlabel("Normalized Frequency  (\times\pi rad/sample)");
ylabel('Power/frequency (W/(rad/sample))');
title("Periodogram (with Hamming window) of centred and detrended signal");
subplot(313);
plot(f_w_original, pxx_w_log_unbiased);
xlabel("Normalized Frequency  (\times\pi rad/sample)");
title("Periodogram (with Hamming window) of centred log signal");
ylabel('Power/frequency (W/(rad/sample))');
