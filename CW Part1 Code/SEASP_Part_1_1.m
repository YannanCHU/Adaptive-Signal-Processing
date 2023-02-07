clc;
clear;
close all;

%% Signal simulation
Fs = 40;            % sampling frequency
Ts = 1 / Fs;        % sampling interval
Ns = 2048;          % number of samples
ts = Ts * (0:1:Ns-1);          % sampling time

% 1) simulate the sinusoidal signal - x1
fx1 = 2;                       % signal frequency
phx1 = 0;                      % phase shift
x1 = 0.8 * cos(2*pi*fx1*ts);   % simulated signal
% figure(); plot(ts, x1); title("The simulated sinusoidal signal");

% 2) simulate the impulse signal - x2
x2 = zeros(1, Ns);
x2(round(Ns/2)) = 1;
% figure(); plot(ts, x2); title("The simulated impulse signal");

%% First definition: Power Spectral Density (PSD) is defined as the DTFT of the ACF
Rxx1 = xcorr(x1, "biased");
Rxx2 = xcorr(x2, "biased");
% index of above calculated Rxx value
ks = -(Ns-1):1:(Ns-1);
% Power Spectral Density (PSD) of two simulated signals
Pxx1_def1 = abs(fftshift(fft(Rxx1)));
Pxx2_def1 = abs(fftshift(fft(Rxx2)));
% The list of frequencies corresponding to the ACF
f_def1 = Fs * ks / (2 * Ns);

%% Second definition: PSD is defined as the expectation of periodogram
% Power Spectral Density (PSD) of two simulated signals
Pxx1_def2 = ( abs(fftshift(fft(x1))) .^ 2 ) / Ns;
Pxx2_def2 = ( abs(fftshift(fft(x2))) .^ 2 ) / Ns;
% The list of frequencies corresponding to the the second definition
f_def2 = linspace(-Fs/2,Fs/2-Fs/(Ns),Ns);

% periodogram function can achieve same results
period1 = periodogram(x1,[],Ns,'centered') * 2 * pi;
period2 = periodogram(x2,[],Ns) * 2 * pi;

%% Comparison between two PSD definitions
figure();
subplot(321);
plot(ts, x1, 'LineWidth', 2);
title("Sinusoidal signal",'FontSize',12);
xlabel("Time (s)",'FontSize',12); ylabel("Signal value",'FontSize',12);
xlim([0, max(ts)]);

subplot(322);
plot(ts, x2, 'LineWidth', 2);
title("Delayed ideal impulse signal",'FontSize',12);
xlabel("Time (s)",'FontSize',12); ylabel("Signal value",'FontSize',12);
xlim([0, max(ts)]);

subplot(323);
plot(ks, Rxx1, 'LineWidth', 2);
title("ACF results of sinusoidal signal (delay slowly)",'FontSize',12);
xlabel("Lag index - k",'FontSize',12); ylabel("ACF- Rxx(k)",'FontSize',12);
xlim([min(ks)-10, max(ks)+10]);

subplot(324);
plot(ks, Rxx2, 'LineWidth', 2);
title("ACF results of impulse signal (delay rapidly)",'FontSize',12);
xlabel("Lag index - k",'FontSize',12); ylabel("ACF- Rxx(k)",'FontSize',12);
xlim([min(ks)-10, max(ks)+10]);

subplot(325);
plot(f_def1, Pxx1_def1, '-', f_def2, Pxx1_def2, '--', 'LineWidth', 2);
title("Two PSD results of sinusoid",'FontSize',12);
legend("Definition 1", "Definition 2",'FontSize',12);
xlabel("Frequency (Hz)",'FontSize',12); ylabel("PSD - P(\omega)",'FontSize',12);
ylim([0, max(Pxx1_def1)+10])

subplot(326);
plot(f_def1, Pxx2_def1, '-', f_def2, Pxx2_def2, '--', 'LineWidth', 2);
title("Two PSD results of impulse",'FontSize',12);
legend("Definition 1", "Definition 2",'FontSize',12);
xlabel("Frequency (Hz)",'FontSize',12); ylabel("PSD - P(\omega)",'FontSize',12);

