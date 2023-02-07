% calculates both biased and unbiased ACF estimates of a signal and then use these 
% ACF estimates to compute the corresponding correlogram

clc;
clear;
close all;

seed = 10;
rng(seed);

%% Signal simulation
Fs = 40;                % sampling frequency
Ts = 1 / Fs;            % sampling interval
Ns = 2048;              % number of samples
ts = Ts * (0:1:Ns-1);   % sampling time

% 1) WGN
WGN = randn(1, Ns);
% 2) Noisy sinusoidal signals
fx1 = 2;    fx2 = 12;                % signal frequency
amp1 = 0.8; amp2 = 1;                % signal amplitude
noisySine = amp1 * cos(2*pi*fx1*ts) + amp2 * cos(2*pi*fx2*ts) + WGN;    % phase = 0
% 3) Filtered WGN
% construct a moving average filter
windowSize = 4; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
filteredWGN = filter(b,a,WGN);

% plot three signals
figure(1);
subplot(211); plot(ts, WGN, ts, filteredWGN, '--'); title("Original and Filtered White Gaussian Noise"); 
xlabel("Time (s)"); ylabel("Signal value"); xlim([0, max(ts)]);
legend("Original", "Filtered");
subplot(212); plot(ts, noisySine); title("Noisy Sinusoidal Signal"); 
xlabel("Time (s)"); ylabel("Signal value"); xlim([0, max(ts)]);

%% ACF calculation - r(k)
% % the custom function - ACF_calulation - can also achieve this goal
% % when 0 <= k <= Ns-1
% % Since Rxx is a symmetric function, the left hand side values can be
% % deduced from the right hand side
% Rxx_WGN_biased2 = ACF_calulation(WGN, "biased");
% Rxx_WGN_unbiased2 = ACF_calulation(WGN, "unbiased");

% 2N-1 ACF values are obtained (-(Ns-1) <= k <= Ns-1)s
Rxx_WGN_biased = xcorr(WGN, "biased");
Rxx_WGN_unbiased = xcorr(WGN, "unbiased");
% compute the correlogram spectral estimator
% Since Pxx is a real function, only the real part of Fourier Transform is
% retained
Pxx_WGN_biased = real(fftshift(fft(ifftshift(Rxx_WGN_biased))));
Pxx_WGN_unbiased = real(fftshift(fft(ifftshift(Rxx_WGN_unbiased))));

Rxx_noisySine_biased = xcorr(noisySine, "biased");
Rxx_noisySine_unbiased = xcorr(noisySine, "unbiased");
% compute the correlogram spectral estimator
Pxx_noisySine_biased = real(fftshift(fft(ifftshift(Rxx_noisySine_biased))));
Pxx_noisySine_unbiased = real(fftshift(fft(ifftshift(Rxx_noisySine_unbiased))));

Rxx_filteredWGN_biased = xcorr(filteredWGN, "biased");
Rxx_filteredWGN_unbiased = xcorr(filteredWGN, "unbiased");
% compute the correlogram spectral estimator
Pxx_filteredWGN_biased = real(fftshift(fft(ifftshift(Rxx_filteredWGN_biased))));
Pxx_filteredWGN_unbiased = real(fftshift(fft(ifftshift(Rxx_filteredWGN_unbiased))));

% plot the ACF results (0 <= k <= Ns-1)
figure(2);
ks = 0:1:Ns-1;
subplot(321);   plot(ks, Rxx_WGN_unbiased(Ns:2*Ns-1), ...
                     ks, Rxx_WGN_biased(Ns:2*Ns-1), 'LineWidth',2);
title("Biased and Unbiased ACF of WGN (0 ≤ k ≤ N-1)",'FontSize',12);
legend("Unbiased", "Biased",'FontSize',11);
xlabel("Lag index - k",'FontSize',12); ylabel("ACF- Rxx(k)",'FontSize',12);
xlim([min(ks), max(ks)]);

subplot(323);   plot(ks, Rxx_noisySine_unbiased(Ns:2*Ns-1), ...
                     ks, Rxx_noisySine_biased(Ns:2*Ns-1), 'LineWidth',2);
title("Biased and Unbiased ACF of Noisy Sinusoidal Signal (0 ≤ k ≤ N-1)",'FontSize',12);
legend("Unbiased", "Biased",'FontSize',11);
xlabel("Lag index - k",'FontSize',12); ylabel("ACF- Rxx(k)",'FontSize',12);
xlim([min(ks), max(ks)]);

subplot(325);   plot(ks, Rxx_filteredWGN_unbiased(Ns:2*Ns-1), ...
                     ks, Rxx_filteredWGN_biased(Ns:2*Ns-1), 'LineWidth',2);
title("Biased and Unbiased ACF of Filtered WGN (0 ≤ k ≤ N-1)",'FontSize',12);
legend("Unbiased", "Biased",'FontSize',11);
xlabel("Lag index - k",'FontSize',12); ylabel("ACF- Rxx(k)",'FontSize',12);
xlim([min(ks), max(ks)]);


% plot the correlogram results
% When the 'unbiased' mode is used, the PSD value may be negative. Thus, the
% log-transfrom cannot be used and PSD cannot be expressed in the unit of dB.
fs = (Fs / Fs) * (-Ns+1:1:Ns-1) / (2 * Ns);
subplot(322);   plot(fs, Pxx_WGN_unbiased, ...
                     fs, Pxx_WGN_biased, 'LineWidth',2);
title("Biased and Unbiased Correlogram of WGN",'FontSize',12);
legend("Unbiased", "Biased",'FontSize',11);
xlabel("Normalized Frequency  (\times\pi rad/sample)",'FontSize',12);
ylabel("Correlogram- Pxx(\omega)",'FontSize',12); 
xlim([min(fs), max(fs)]);

subplot(324);   plot(fs, Pxx_noisySine_unbiased, ...
                     fs, Pxx_noisySine_biased, 'LineWidth',2);
title("Biased and Unbiased Correlogram of Noisy Sinusoidal Signal",'FontSize',12);
legend("Unbiased", "Biased",'FontSize',11);
xlabel("Normalized Frequency  (\times\pi rad/sample)",'FontSize',12); 
ylabel("Correlogram- Pxx(\omega)",'FontSize',12);
xlim([min(fs), max(fs)]);

subplot(326);   plot(fs, Pxx_filteredWGN_unbiased, ...
                     fs, Pxx_filteredWGN_biased, 'LineWidth',2);
title("Biased and Unbiased Correlogram of Filtered WGN",'FontSize',12);
legend("Unbiased", "Biased",'FontSize',11);
xlabel("Normalized Frequency  (\times\pi rad/sample)",'FontSize',12);
ylabel("Correlogram- Pxx(\omega)",'FontSize',12);
xlim([min(fs), max(fs)]);


%%

function Rxx = ACF_calulation(x, flag)
    % Biased ACF is used by default
    % x: original signal samples
    % flag: either 'biased' or 'unbiased'
    N = length(x);
    Rxx = zeros(1,N);

    if strcmp(flag, 'unbiased')
        for k=0:N-1
            % normalisation term for unbiased ACF
            gain = 1/(N-k);
            Rxx(k+1) = gain * sum(x(k+1:N).*conj(x(1:N-k)));
        end
    else
        % normalisation term for biased ACF
        gain = 1/N;
        for k=0:N-1
            Rxx(k+1) = gain * sum(x(k+1:N).*conj(x(1:N-k)));
        end
    end
end