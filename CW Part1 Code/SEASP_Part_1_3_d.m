% MUSIC Algorithm
% Complex sinusoidal signal + complex noise generation
clc;
clear;
close all;

seed = 3;
rng(seed);

%% Signal simulation and PSD estimation
Fs = 1;                     % sampling frequency (unit sampling rate)
Ts = 1 / Fs;                % sampling interval
Ns = [30, 45, 70, 120];     % number of samples
nffts = 4*128;              % number of DFT points

% store the PSD estimation results
Pxxs = zeros(length(Ns), nffts);

for i = 1:length(Ns)
    % sampling time
    ts = Ts * (0:1:Ns(i)-1);
    % complex noise
    noiseVar = 0.2;                       % noise variance
    noise = noiseVar/sqrt(2)*(randn(size(ts))+1j*randn(size(ts)));
    % complex noisy signal
    fx1 = 0.3;    fx2 = 0.32;             % signal frequency
    amp1 = 1;     amp2 = 1;               % signal amplitude
    x = amp1*exp(1j*2*pi*fx1*ts)+ ...
        amp2*exp(1j*2*pi*fx2*ts)+ noise;

    [Pxxs(i,:),f] = periodogram(x,rectwin(length(x)),nffts,Fs);
end


figure();
plot(1000*f, pow2db(Pxxs),'LineWidth', 2);
xlim([0, 700]);
title('Periodogram (with rectangular window) of noisy complex exponential signal','FontSize',14);
xlabel('Frequency (mHz)','FontSize',14);
ylabel('Power/frequency (dB/Hz)','FontSize',14);
legend(sprintf("N = %d", Ns(1)), sprintf("N = %d", Ns(2)), ...
    sprintf("N = %d", Ns(3)), sprintf("N = %d", Ns(4)),'FontSize',12);
grid on;
