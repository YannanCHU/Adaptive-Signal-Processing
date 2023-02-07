% MUSIC Algorithm
% Complex sinusoidal signal + complex noise generation
clc;
clear;
close all;

seed = 12;
rng(seed);

%% Signal simulation and PSD estimation
Fs = 1;                     % sampling frequency (unit sampling rate)
Ts = 1 / Fs;                % sampling interval
Ns = 30;                    % number of samples
numOfRealisations = 500;    % The number of realisations
nffts = 2*128;              % number of DFT points

% create a matrix to store Pseudospectrum estimate of each trial
PseudoSpecs = zeros(numOfRealisations, nffts);

for i = 1:numOfRealisations
    % sampling time
    ts = Ts * (0:1:Ns-1);
    % complex noise
    noiseVar = 0.2;                       % noise variance
    noise = noiseVar/sqrt(2)*(randn(size(ts))+1j*randn(size(ts)));
    % complex noisy signal
    fx1 = 0.3;    fx2 = 0.32;             % signal frequency
    amp1 = 1;     amp2 = 1;               % signal amplitude
    x = amp1*exp(1j*2*pi*fx1*ts)+ ...
        amp2*exp(1j*2*pi*fx2*ts)+ noise;

    [X,R] = corrmtx(x,14,'modified');
    [PseudoSpecs(i,:),F] = pmusic(R,2,[ ],1,'corr');
end

figure(1);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
p1(1:numOfRealisations) = plot(F,PseudoSpecs, 'c','linewidth',1); set(gca,'xlim',[0.25 0.40]);
grid on; xlabel('Frequency (Hz)', 'FontSize',14);  ylabel('Pseudospectrum', 'FontSize',14);
title("Mean pseudospectrum of "+numOfRealisations+" trials" + ...
    " (N = "+Ns+")", 'FontSize',14); hold on;
p1(numOfRealisations+1) = plot(F,mean(PseudoSpecs), 'b','linewidth',2); hold off;
legend(p1([1,numOfRealisations+1]), 'realisations', 'mean');

nexttile;
plot(F,std(PseudoSpecs), 'r', 'linewidth',2); set(gca,'xlim',[0.25 0.40]);
grid on; xlabel('Frequency (Hz)', 'FontSize',14);  ylabel('Pseudospectrum', 'FontSize',14);
title("Standard deviation of pseudospectrum of "+numOfRealisations+" trials" + ...
    " (N = "+Ns+")", 'FontSize',14);

%%  noiseVar = 0.3;
noiseVarH = 0.3;     % noise variance
PseudoSpecsHVar = zeros(numOfRealisations, nffts);
for i = 1:numOfRealisations
    % sampling time
    ts = Ts * (0:1:Ns-1);
    % complex noise
    
    noise = noiseVarH/sqrt(2)*(randn(size(ts))+1j*randn(size(ts)));
    % complex noisy signal
    fx1 = 0.3;    fx2 = 0.32;             % signal frequency
    amp1 = 1;     amp2 = 1;               % signal amplitude
    x = amp1*exp(1j*2*pi*fx1*ts)+ ...
        amp2*exp(1j*2*pi*fx2*ts)+ noise;

    [X,R] = corrmtx(x,14,'modified');
    [PseudoSpecsHVar(i,:),F] = pmusic(R,2,[ ],1,'corr');
end

%%  noiseVar = 0.1;
noiseVarL = 0.1;     % noise variance
PseudoSpecsLVar = zeros(numOfRealisations, nffts);
for i = 1:numOfRealisations
    % sampling time
    ts = Ts * (0:1:Ns-1);
    % complex noise
    
    noise = noiseVarL/sqrt(2)*(randn(size(ts))+1j*randn(size(ts)));
    % complex noisy signal
    fx1 = 0.3;    fx2 = 0.32;             % signal frequency
    amp1 = 1;     amp2 = 1;               % signal amplitude
    x = amp1*exp(1j*2*pi*fx1*ts)+ ...
        amp2*exp(1j*2*pi*fx2*ts)+ noise;

    [X,R] = corrmtx(x,14,'modified');
    [PseudoSpecsLVar(i,:),F] = pmusic(R,2,[ ],1,'corr');
end

%% Analyse the effect of noise level

figure();
subplot(211);
plot(F,mean(PseudoSpecsHVar), ...
    F,mean(PseudoSpecs), ...
    F,mean(PseudoSpecsLVar), ...
    'linewidth',2); 
set(gca,'xlim',[0.25 0.40]);
grid on; xlabel('Frequency (Hz)');  ylabel('Pseudospectrum');
title("Mean pseudospectrum of "+numOfRealisations+" trials" + ...
    " (N = "+Ns+")"); hold on;
legend("Noise Variance = " + noiseVarH,...
       "Noise Variance = " + noiseVar,...
       "Noise Variance = " + noiseVarL);

subplot(212);
plot(F,std(PseudoSpecsHVar), ...
    F,std(PseudoSpecs), ...
    F,std(PseudoSpecsLVar), ...
    'linewidth',2); set(gca,'xlim',[0.25 0.40]);
grid on; xlabel('Frequency (Hz)');  ylabel('Pseudospectrum');
title("Standard deviation of pseudospectrum of "+numOfRealisations+" trials" + ...
    " (N = "+Ns+")");
legend("Noise Variance = " + noiseVarH,...
       "Noise Variance = " + noiseVar,...
       "Noise Variance = " + noiseVarL);

