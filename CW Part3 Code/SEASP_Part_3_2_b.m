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
figure(1); subplot(2,1,1); plot(n, fn, 'LineWidth', 2); 
title("Frequency of FM signal - f(n)", 'FontSize', 14); 
xlabel("Sample index - n", 'FontSize', 14); ylabel("f(n) (Hz)", 'FontSize', 14);

subplot(2,1,2); plot(n, angle(exp(1j*( 2*pi*Phi(n) / fs ))), 'LineWidth', 2); 
title("Phase of FM signal - \phi(n)", 'FontSize', 14)
xlabel("Sample index - n", 'FontSize', 14); ylabel("\phi(n) (rad)", 'FontSize', 14);

rng(12);
eta = sqrt(noiseVar / 2) * (randn(1,N) + 1j * randn(1,N));    % the noise signal
coefficient = circularityMeasure(eta);
disp("The circularity of complex noise is " + coefficient);
yn = exp(1j*( 2*pi*Phi(n) / fs )) + eta;    % the FM signal

%%
mus = [0.003, 0.01, 0.03, 0.3];      % step size
numOfPoints = 1024;
CLMS_sepctrums = cell(length(mus),1);

for k = 1:1:length(mus)
    H = zeros(numOfPoints, N);
    [~, ~, AR_coefs] = CLMS_arma([], yn, mus(k), AR_order, 0);
    for i = 1:N
        %  Run complex-valued LMS algorithm to estimate AR coefficient
        [hFreResp, wAngFreResp] = freqz(1, [1; -conj(AR_coefs(i))], numOfPoints, fs);  % compute power spectrum
        H(:,i) = abs(hFreResp) .^ 2;    % store it in a matrix
    end

    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    CLMS_sepctrums{k,1} = H;
end

figure(2);
tiledlayout(2,2,'TileSpacing','compact');

nexttile;
surf(n, wAngFreResp, CLMS_sepctrums{1,1},'LineStyle','none');     view(2);
title("CLMS based spectrum estimate (\mu = " + mus(1) +")", 'FontSize',14);    c = colorbar;   ylabel(c, "PSD", 'FontSize',14);
xlabel("Sample Index", 'FontSize',14);
ylabel("Frequency (Hz)", 'FontSize',14);
ylim([0, max(wAngFreResp)]);

nexttile;
surf(n, wAngFreResp, CLMS_sepctrums{2,1},'LineStyle','none');     view(2);
title("CLMS based spectrum estimate (\mu = " + mus(2) +")", 'FontSize',14);    c = colorbar;   ylabel(c, "PSD", 'FontSize',14);
xlabel("Sample Index", 'FontSize',14);
ylabel("Frequency (Hz)", 'FontSize',14);
ylim([0, max(wAngFreResp)]);

nexttile;
surf(n, wAngFreResp, CLMS_sepctrums{3,1},'LineStyle','none');     view(2);
title("CLMS based spectrum estimate (\mu = " + mus(3) +")", 'FontSize',14);    c = colorbar;   ylabel(c, "PSD", 'FontSize',14);
xlabel("Sample Index", 'FontSize',14);
ylabel("Frequency (Hz)", 'FontSize',14);
ylim([0, max(wAngFreResp)]);

nexttile;
surf(n, wAngFreResp, CLMS_sepctrums{4,1},'LineStyle','none');     view(2);
title("CLMS based spectrum estimate (\mu = " + mus(4) +")", 'FontSize',14);    c = colorbar;   ylabel(c, "PSD", 'FontSize',12);
xlabel("Sample Index", 'FontSize',14);
ylabel("Frequency (Hz)", 'FontSize',14);
ylim([0, max(wAngFreResp)]);

