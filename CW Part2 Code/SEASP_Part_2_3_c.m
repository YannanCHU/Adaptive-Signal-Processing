clc;
clear;
close all;

%% signal generation
% clean signal generation
w0 = 0.01 * pi;     % angular frequency
N = 1000;           % number of samples
t = 0:1:N-1;
x = (sin(w0*t)).';

% colored noise generation
rng(0);
noiseVar = 1;
numOfRealisations = 100;
a1 = [0, 0.5];
MA_MDL = arima('MA', a1, 'Variance', noiseVar, 'Constant', 0);
% v is the white noise while eta is the colored noise
[eta, v] = simulate(MA_MDL, N, 'NumPaths', numOfRealisations);

%  noise-corrupted signal generation
s = x + eta;

% reference noise - eta1
eta1 = 0.8 * eta + 0.1;

%% Calculate the prediction error
delay = 3;                  % delay
M = 5;                      % order of filter
mu = 0.01;                  % step size
steadyStartIndex = 300;     % the index within the steady state range

MSPE_ALE = 0;
MSPE_ANC = 0;
MSPEs_ALE_temp = zeros(numOfRealisations,1);
MSPEs_ANC_temp = zeros(numOfRealisations,1);
x_hat_ALE = cell(numOfRealisations, 1);
x_hat_ANC = cell(numOfRealisations, 1);

for j = 1:numOfRealisations
    % the ALE results
    [x_hat_ALE{j,1}, ~, ~] = LMS_ale_anc([], s(:,j), mu, M, delay, 'ALE');
    MSPEs_ALE_temp(j) = mean(( x(steadyStartIndex:end) - x_hat_ALE{j,1}(steadyStartIndex:end) ).^2);
    
    % the ANC results
    [noiseEst, ~, ~] = LMS_ale_anc(eta1(:,j), s(:,j), mu, M, [], 'ANC');
    x_hat_ANC{j,1} = s(:,j) - noiseEst;
    MSPEs_ANC_temp(j) = mean(( x(steadyStartIndex:end) - x_hat_ANC{j,1}(steadyStartIndex:end) ).^2);
end

MSPEs_ALE = pow2db(mean(MSPEs_ALE_temp));
MSPEs_ANC = pow2db(mean(MSPEs_ANC_temp));

%% display the results
figure(1);
tiledlayout(3,1,'TileSpacing','compact');
nexttile;
p1 = plot(t, s, 'c');   hold on;
p2 = plot(t, cell2mat(x_hat_ALE.'), 'g');
p3 = plot(t, x, 'r', 'LineWidth', 2);   hold off;
legend([p1(1), p2(1), p3], ...
    {'noisy $s$', 'ALE $\hat{x}$','clean $x$'}, ...
    'Interpreter','latex', 'FontSize',12);
title("ALE denoising results (\Delta = " + delay + ", M = " + M + ", MSPE = "+num2str(MSPEs_ALE,3)+" dB)", 'FontSize',14);
xlabel("Sample Index", 'FontSize',14);
ylabel("Amplitude", 'FontSize',14);

nexttile;
p1 = plot(t, s, 'c');   hold on;
p2 = plot(t, cell2mat(x_hat_ANC.'), 'color', '#EDB120');
p3 = plot(t, x, 'r', 'LineWidth', 2);   hold off;
legend([p1(1), p2(1), p3], ...
    {'noisy $s$', 'ANC $\hat{x}$','clean $x$'}, ...
    'Interpreter','latex', 'FontSize',12);
title("ANC denoising results (M = " + M + ", MSPE = "+num2str(MSPEs_ANC,3)+" dB)", 'FontSize',14);
xlabel("Sample Index", 'FontSize',14);
ylabel("Amplitude", 'FontSize',14);

% figure(2);
nexttile;
p1 = plot(t, mean(s,2), 'c', 'LineWidth', 2);   hold on;
p2 = plot(t, mean(cell2mat(x_hat_ALE.'), 2), 'g', 'LineWidth', 2);
p3 = plot(t, mean(cell2mat(x_hat_ANC.'), 2), 'color', '#EDB120', 'LineWidth', 2);
p4 = plot(t, x, 'r', 'LineWidth',2);    hold off;
title("The ALE and ANC denoising results", 'FontSize',14)
xlabel("Sample Index", 'FontSize',14);
ylabel("Amplitude", 'FontSize',14);
legend([p1(1), p2(1), p3(1), p4], ...
    {'noisy $s$', 'ALE $\hat{x}$', 'ANC $\hat{x}$','clean $x$'}, ...
    'Interpreter','latex', 'FontSize',12);


