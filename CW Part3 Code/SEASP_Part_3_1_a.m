clc;
clear;
close all;

% coefficients of WLMA
b1 = 1.5 + 1i;
b2 = 2.5 - 0.5i;
b_coef = [b1, b2];
MA_order = length(b_coef);

% parameter of the noise signal
N = 1000;
numOfRealisations = 100;

% create the complex white Gaussian noise
rng(0);
noiseVar = 1;
x = (noiseVar / sqrt(2)) * (randn(numOfRealisations, N) + 1i * randn(numOfRealisations, N));
x1 = [zeros(numOfRealisations, 1), x(:, 1:end-1)];

% create the moving averaged complex signal
y = x + b1 * x1 + b2 * conj(x1);

% Circularity coefficient
coef_x = circularityMeasure(x(:).');
coef_y = circularityMeasure(y(:).');

%% apply CLMS to identify the WLMA model
mu = 0.08;      % step size
errors_CLMS = cell(1,numOfRealisations);
errors_ACLMS = cell(1,numOfRealisations);
for i = 1:1:numOfRealisations
    [~, error_clms, ~] = CLMS_arma(x(i,:), y(i,:), mu, 0, MA_order);
    errors_CLMS{1,i} = (abs(error_clms).^2);
    [~, error_aclms, ~, ~] = ACLMS_arma(x(i,:), y(i,:), mu, 0, MA_order);
    errors_ACLMS{1,i} = (abs(error_aclms).^2);
end

% mean(cell2mat(errors_CLMS),2);
% mean(cell2mat(errors_ACLMS),2);


%% display results
figure(1);
% tiledlayout(2,1,'TileSpacing','compact');
% nexttile;
scatter(real(y(:)), imag(y(:)), 96, '.'); hold on;
scatter(real(x(:)), imag(x(:)), 96, '.'); hold off;
daspect([1,1,1]);
axis([-13 13 -6 6]);
title("Real and imaginary parts of WGN and WLMA(1) process",'FontSize',16);
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("WLMA(1) (|\rho| = " +coef_y+")", ...
    "WGN (|\rho| = " +coef_x+")", 'FontSize',14,'Location','southeast');

figure(2);
% nexttile;
plot(1:N, pow2db(mean(cat(3,errors_CLMS{:}),3)), 1:N, pow2db(mean(cat(3,errors_ACLMS{:}),3)), 'LineWidth', 2);
title("Learning curves for ACLMS and CLMS (Step size \mu = "+mu+")", 'FontSize',16);
xlabel("Sample Index", 'FontSize',14);
ylabel("Squared prediction error in dB", 'FontSize',14);
legend("CLMS", "ACLMS", 'FontSize', 12, 'Location','east');
