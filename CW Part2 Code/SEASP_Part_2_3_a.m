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
% rng(0);
% noiseVar = 1;
% v = sqrt(noiseVar) * randn(1,N);     % white noise
% v_2_delay = [0, 0, v(1:end-2)];
% eta = v + 0.5 * v_2_delay;           % colored noise

rng(0);
noiseVar = 1;
numOfRealisations = 100;
a1 = [0, 0.5];
MA_MDL = arima('MA', a1, 'Variance', noiseVar, 'Constant', 0);
% v is the white noise while eta is the colored noise
[eta, v] = simulate(MA_MDL, N, 'NumPaths', numOfRealisations);

%  noise-corrupted signal generation
s = x + eta;

%% Calculate the prediction error
delays = 1:1:10;
filterOrder = 5;
mu = 0.01;
MSPE = zeros(length(delays),1);
MSPEs_temp = zeros(length(delays), numOfRealisations);
x_hat = cell(length(delays), numOfRealisations);
steadyStartIndex = 300;

for i = 1:1:length(delays)
    delay = delays(i);
    
    for j = 1:numOfRealisations
        [x_hat{i,j}, ~, ~] = LMS_ale(s(:,j), mu, filterOrder, delay);
        MSPEs_temp(i,j) = mean(( x(steadyStartIndex:end) - x_hat{i,j}(steadyStartIndex:end) ).^2);
    end

    MSPE(i) = mean(MSPEs_temp(i,:), 2);
end

%% results demonstration

figure(1);
plot(delays, pow2db(MSPE), 'LineWidth', 2); title(["Mean Square Prediction Error against delay" + " (M = " + filterOrder+")"]);
xlabel("Delay"); ylabel("MSPE (dB)");

figure(2);
tiledlayout(1,4,'TileSpacing','compact');
for delay_idx = 1:1:4
    nexttile;
    p1 = plot(t, s, 'b');   hold on;
    p2 = plot(t, cell2mat(x_hat(delay_idx,:)), 'g');
    p3 = plot(t, x, 'r', 'LineWidth',2);    hold off;
    title(["ALE results when M = "+filterOrder+", \Delta = " + delays(delay_idx), " (MSPE = "+num2str(pow2db(MSPE(delay_idx)),3)+" dB)"], 'FontSize',13);
    xlabel("Smaple Index", 'FontSize',13); ylabel("Amplitude", 'FontSize',13);  
end
legend([p1(1), p2(1), p3], {'noisy $s$', 'denoised $\hat{x}$', 'clean $x$'}, 'Interpreter','latex', 'FontSize',13);
