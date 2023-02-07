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
delays = 3:1:25;    % different delays
filterOrders = 5:5:20;    % order of filter
mu = 0.01;      % step size
MSPE = zeros(length(delays), length(filterOrders));
MSPEs_temp = zeros(length(delays), numOfRealisations, length(filterOrders));
x_hat = cell(length(delays), numOfRealisations, length(filterOrders));
steadyStartIndex = 300;     % the index within the steady state range

for k = 1:1:length(filterOrders)
    M = filterOrders(k);
    for i = 1:1:length(delays)
        delay = delays(i);

        for j = 1:numOfRealisations
            [x_hat{i,j,k}, ~, ~] = LMS_ale(s(:,j), mu, M, delay);
            MSPEs_temp(i,j,k) = mean(( x(steadyStartIndex:end) - x_hat{i,j,k}(steadyStartIndex:end) ).^2);
        end
        
        MSPE(i,k) = mean(MSPEs_temp(i,:,k), 2);
    end
end
%% results demonstration

% MSPE vrs. delay
figure(1);
subplot(2,1,1);
plot(delays, pow2db(MSPE.'), 'LineWidth', 2); title("Mean Square Prediction Error against delay", 'FontSize',16);
xlabel("Delay", 'FontSize',14); ylabel("MSPE (dB)", 'FontSize',14);
legend("M = " + filterOrders(1), "M = " + filterOrders(2), "M = " + filterOrders(3), "M = " + filterOrders(4), 'Location','southeast', 'FontSize',12);
xlim([min(delays), max(delays)]);
%% Mean of estimated signal vrs. smaple index
figure(2);
tiledlayout(1,4,'TileSpacing','compact');
order_idx = 1;
for delay_idx = 1:7:25
    nexttile;
    p1 = plot(t, mean(s,2), 'b');   hold on;
    p2 = plot(t, mean(cell2mat(x_hat(delay_idx,:,1)),2), 'g');
    p3 = plot(t, x, 'r', 'LineWidth',2);    hold off;
    title(["The ALE results when M = "+filterOrders(order_idx)+", \Delta = " + delays(delay_idx), " (MSPE = "+num2str(pow2db(MSPE(delay_idx, order_idx)),3)+" dB)"]);
    xlabel("Smaple Index"); ylabel("Amplitude");  
end
legend([p1(1), p2(1), p3], {'noisy $s$', 'denoised $\hat{x}$', 'clean $x$'}, 'Interpreter','latex', 'FontSize',12);
%% MSPE vrs. order of filter
delays2 = [3,15,25];
filterOrders2 = 1:1:25;
MSPE2 = zeros(length(delays2), length(filterOrders2));
MSPEs_temp2 = zeros(length(delays2), numOfRealisations, length(filterOrders2));
x_hat2 = cell(length(delays2), numOfRealisations, length(filterOrders2));
for k = 1:1:length(filterOrders2)
    M = filterOrders2(k);
    for i = 1:1:length(delays2)
        delay = delays2(i);

        for j = 1:numOfRealisations
            [x_hat2{i,j,k}, ~, ~] = LMS_ale(s(:,j), mu, M, delay);
            MSPEs_temp2(i,j,k) = mean(( x(steadyStartIndex:end) - x_hat2{i,j,k}(steadyStartIndex:end) ).^2);
        end
        
        MSPE2(i,k) = mean(MSPEs_temp2(i,:,k), 2);
    end
end

figure(1);
subplot(2,1,2);
plot(filterOrders2, pow2db(MSPE2), 'LineWidth',2); title("Mean Square Prediction Error against order of filter", 'FontSize',16);
xlabel("Order of filter", 'FontSize',14); ylabel("MSPE (dB)", 'FontSize',14);
legend("\Delta = "+delays2(1), "\Delta = "+delays2(2), "\Delta = "+delays2(3), 'FontSize',12, 'Location','southeast');