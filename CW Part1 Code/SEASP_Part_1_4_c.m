% Auto-regression model
% Analysis on the order of AR model
clc;
clear;
close all;

seed = 0;
rng(seed);

%% Signal simulation
Fs = 1;                     % sampling frequency (unit sampling rate)
Ts = 1 / Fs;                % sampling interval
Ns = 10000;                 % number of samples
coef = [2.76, -3.81, 2.65, -0.92];  % coefficients of AR model
noiseVar = 1;               % noise variance

% creat an array to store Ns elements
x = zeros(Ns,1);
for n = 1+length(coef):Ns
    x(n) = coef(1)*x(n-1) + coef(2)*x(n-2) + coef(3)*x(n-3) + coef(4)*x(n-4) + noiseVar*randn(1,1);
end

% discard the first 500 samples
x = x(501:end);

%% true PSD of AR model
[ht,f] = freqz(1,[1, -coef],Ns,Fs);   % get the frequency reeponse of filter 
PxxTrue = pow2db(abs(ht).^2);        % get the true PSD

%% estimated PSDs of different orders
% order of underlying AR model
p = 2:1:14;
% create a cell to store estimated AR parameters
ARparameters = cell(length(p),1);
% create a vector to store estimated noise variance
noiseVarss = zeros(length(p),1);
% create a matrix to store estimated PSD
PxxEst = zeros(length(f),length(p));
% create a matrix to store the difference between estimated and true PSD
PxxDiff = zeros(length(f),length(p));

for i = 1:1:length(p)
    % normalized autoregressive (AR) parameters corresponding to a model of
    % order p for the input array x;
    [ARparameters{i},noiseVars(i),~] = aryule(x,p(i));
    % estimated PSD of AR model
    he = freqz(sqrt(noiseVars(i)),ARparameters{i},Ns,Fs);   % get the frequency reeponse of filter
    PxxEst(:,i) = pow2db(abs(he).^2);        % get the true PSD

    % % Autoregressive power spectral density estimate — Yule-Walker method
    % [PxxEst2,f] = pyulear(x,p,Ns*2-1,Fs);

    % difference between estimated and ground truth PSD
    PxxDiff(:,i) = abs(PxxTrue-PxxEst(:,i));
end

%% plot the true and estimated PSD as well as the MSE of each order
figure();
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(f,PxxTrue, f,PxxEst(:,1), f,PxxEst(:,3), f,PxxEst(:,5), 'LineWidth', 2);
title("The true PSD and Estimated PSD against frequency ("+length(x) +" samples are used)", 'FontSize',14);
xlabel("Frequency (Hz)", 'FontSize',14); ylabel("PSD", 'FontSize',14);
legend("True PSD", "order = " + p(1) + " (Low)", "order = " + p(3) + " (Correct)", ...
    "order = " + p(5) + " (High)", 'FontSize', 12)

nexttile;
semilogy(p, mean(PxxDiff.^2), 'LineWidth',2);
title("The MSE between estimated and true PSDs", 'FontSize',14);
xlabel("Order", 'FontSize',14); ylabel("MSE", 'FontSize',14);
