clc;
clear;
close all;

%% generate the frequency modulated (FM) signal - yn
noiseVar = 0.05;    % noise variance
N = 1500;           % number of signals
n = 1:1:N;          % the sequence of smaple indexes
fs = 1500;          % the sampling frequency
AR_order = 1;       % the order of AR filter

fn = (1<=n) .* (n<=500) .* 100 + ...
     (501<=n) .* (n<=1000) .* (100 + (n-500)/2) + ...
     (1001<=n) .* (n<=1500) .* (100 + ((n-1000)/25).^2);

Phi = cumsum(fn);   % the phase signal
figure(1); plot(n, fn, 'LineWidth', 2); 
title("Frequency of FM signal - f(n)", 'FontSize', 14); 
xlabel("Sample index - n", 'FontSize', 14); ylabel("f(n) (Hz)", 'FontSize', 14);
ylim([90, 510]);

figure(2); plot(n, angle(exp(1j*( 2*pi*Phi(n) / fs ))), 'LineWidth', 2); 
title("Phase of FM signal - \phi(n)", 'FontSize', 14)
xlabel("Sample index - n", 'FontSize', 14); ylabel("\phi(n) (rad)", 'FontSize', 14);

rng(12);
eta = sqrt(noiseVar / 2) * (randn(1,N) + 1j * randn(1,N));    % the noise signal
coefficient = circularityMeasure(eta);
disp("The circularity of complex noise is " + coefficient);
yn = exp(1j*( 2*pi*Phi(n) / fs )) + eta;    % the FM signal

%% plot the spectrum and analyse the effect of AR filter order
AR_coefs = aryule(yn, AR_order);
[hFreResp, wAngFreResp] = freqz(1, AR_coefs, N, fs);

figure(3); 
lineWidths = 2*ones(10,1); lineWidths(1) = 4;
for order_ith = 1:3:10
    AR_coefs = aryule(yn, order_ith);
    [hFreResp, wAngFreResp] = freqz(1, AR_coefs, N, fs);
    plot(wAngFreResp, pow2db(abs(hFreResp).^2), 'LineWidth', lineWidths(order_ith)); hold on;
end
hold off;
title("Power spectrum of the whole FM signal with different AR filter orders", 'FontSize', 14);
xlabel("Frequency (Hz)", 'FontSize', 14);
ylabel("Magnitude (dB)", 'FontSize', 14);
legend("AR("+[1:3:10]+")", 'FontSize', 12);
xlim([0,max(wAngFreResp)]);
grid on;
