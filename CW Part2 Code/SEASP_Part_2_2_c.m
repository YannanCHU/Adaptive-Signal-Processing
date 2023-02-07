% three GASS (gradient adaptive step-size) algorithms

clc;
clear;
close all;

% number of samples
N = 1000;
% step sizes
mus = [0.1, 2];
% coefficients of MA model
a1 = [0.9];
% order of MA model
MA_order = length(a1);
% noise variance
noiseVar = 0.5;
% number of realisations
numOfRealisations = 100;

% get the AR(1) model
% x(n) = a1*noise(n-1) + noise(n)
% rng(0);
% noise = sqrt(noiseVar) * randn(N,1);
% noise_1_delay = [0; noise(1:end-1)];
% x = 0.9*noise_1_delay + noise;

rng(0);
MA_MDL = arima('MA', a1, 'Variance', noiseVar, 'Constant', 0);
[x100, noise1] = simulate(MA_MDL, N, 'NumPaths', numOfRealisations);

%% GNGD algorithm
% store the weight error - Wo - W(n)
weighErrorGNGD = cell(length(mus), numOfRealisations);
averageWeightErrorGNGD = cell(length(mus),1);
weightsGNGD = cell(length(mus), numOfRealisations);

errorGNGD_sq = cell(length(mus), numOfRealisations);
averageErrorGNGD = cell(length(mus),1);
averageWeightsGNGD = cell(length(mus),1);

rho_gngd = 5e-3;

for i = 1:length(mus)
    for j = 1:numOfRealisations
        % the input signal is the noise
        % thus, the real order of MA should + 1
        [~, errorLMS, weights] = GNGD_arma(noise1(:,j), x100(:,j), mus(i), rho_gngd, 0, MA_order+1);
        errorGNGD_sq{i,j} = (errorLMS .^ 2);
        weighErrorGNGD{i,j} = a1 - weights(2,:);
        weightsGNGD{i,j} = weights(2,:);
    end
    averageErrorGNGD{i,1} = 10*log10(mean(cat(3, errorGNGD_sq{i,:}),3));
    averageWeightErrorGNGD{i,1} = mean(cat(3, weighErrorGNGD{i,:}),3);
    averageWeightsGNGD{i,1} = mean(cat(3, weightsGNGD{i,:}),3);
end

%% GASS - Benveniste's algorithm
mu_inits = [0.1 0.54];  % initial guess of the step-size
rho_gass = rho_gngd;     % learning rate for adaptive step-size algorithm

averageErrors_ben = cell(length(mu_inits),1);
averageWeightsError_ben = cell(length(mu_inits),1);
averageWeightsGASS = cell(length(mu_inits),1);

errors_ben = cell(1,numOfRealisations);
weightsError_ben = cell(1,numOfRealisations);
weightsGASS = cell(length(mus), numOfRealisations);

for i = 1:length(mu_inits)
    mu_init = mu_inits(i);
    for j = 1:numOfRealisations
        [~, errors, weights] = GASS_arma(noise1(:,j), x100(:,j), rho_gass, mu_init, 0, MA_order+1, 'Benveniste', []);
        weightsError_ben{j} = (a1 - weights(2,:)).';
        errors_ben{j} = (errors .^ 2);
        weightsGASS{i,j} = weights(2,:);
    end

    averageErrors_ben{i,1} = pow2db(mean(cell2mat(errors_ben),2));
    averageWeightsError_ben{i,1} = mean(cell2mat(weightsError_ben),2);
    averageWeightsGASS{i,1} = mean(cat(3, weightsGASS{i,:}),3);
end

figure();
% subplot(2,1,1);
plot(1:N, a1*ones(1,N), 'g', ...
    1:N, averageWeightsGNGD{1,1}, 'b-', 1:N, averageWeightsGASS{1,1}, 'c--', ...
    1:N, averageWeightsGNGD{2,1}, 'm-', 1:N, averageWeightsGASS{2,1}, 'r--', ...
    'LineWidth', 2);
xlabel("Sample index", 'FontSize',14);
ylabel("Weight Estimates", 'FontSize', 14);
title("Weight estimates of Benveniste’s GASS algorithm and GNGD", 'FontSize',16);
legend("Ground Truth Weight Coefficient", ...
    "GNGD (\mu = "+ mus(1) +", \rho = "+rho_gngd+")", ...
    "Benveniste's GASS (\mu_{initial} = "+mu_inits(1)+", \rho = "+rho_gass+")", ...
    "GNGD (\mu = "+ mus(2) +", \rho = "+rho_gngd+")", ...
    "Benveniste's GASS (\mu_{initial} = "+mu_inits(2)+", \rho = "+rho_gass+")", ...
    'FontSize',10);
xlim([0,200]);
% ylim([-0.1, 1]);

figure();
% subplot(2,1,2);
plot(1:N, averageWeightErrorGNGD{1,1}, 'b-', 1:N, averageWeightsError_ben{1,1}, 'c-', ...
    1:N, averageWeightErrorGNGD{2,1}, 'm-', 1:N, averageWeightsError_ben{2,1}, 'r-','LineWidth', 2);
xlabel("Sample index", 'FontSize',14);
ylabel("Weight Error", 'FontSize', 14);
title("The weight error curve of Benveniste’s GASS algorithm and GNGD", 'FontSize',16);
legend("GNGD (\mu = "+ mus(1) +", \rho = "+rho_gngd+")", "Benveniste's GASS (\mu_{initial} = "+mu_inits(1)+", \rho = "+rho_gass+")", ...
    "GNGD (\mu = "+ mus(2) +", \rho = "+rho_gngd+")", ...
    "Benveniste's GASS (\mu_{initial} = "+mu_inits(2)+", \rho = "+rho_gass+")", 'FontSize',10);
xlim([0,200]);
% ylim([-2,1]);

figure();
% subplot(2,1,2);
plot(1:N, averageErrorGNGD{1,1}, 'b-', 1:N, averageErrors_ben{1,1}, 'c-', ...
    1:N, averageErrorGNGD{2,1}, 'm-', 1:N, averageErrors_ben{2,1}, 'r-','LineWidth', 2);
xlabel("Sample index", 'FontSize',14);
ylabel("Weight Error (dB)", 'FontSize', 14);
title("The squared prediction error curve of Benveniste’s GASS algorithm and GNGD", 'FontSize',16);
legend("GNGD (\mu = "+ mus(1) +", \rho = "+rho_gngd+")", "Benveniste's GASS (\mu_{initial} = "+mu_inits(1)+", \rho = "+rho_gass+")", ...
    "GNGD (\mu = "+ mus(2) +", \rho = "+rho_gngd+")", ...
    "Benveniste's GASS (\mu_{initial} = "+mu_inits(2)+", \rho = "+rho_gass+")", 'FontSize',10);
xlim([0,200]);
% ylim([-350, 350]);

%     figure(2);
%     subplot(length(mu_inits),1,k);
%     plot(1:N, averageErrorGNGD{1,1}, 1:N, averageErrorGNGD{2,1}, ...
%         1:N, averageErrors_ben, 1:N, averageErrors_ang, ...
%         1:N, averageErrors_mat, ...
%         'LineWidth', 2);
%     xlabel("Sample index", 'FontSize',14);
%     ylabel("Squared prediction error in dB", 'FontSize', 14);
%     title("The aveage squared error curve of 100 realisations of x(n)", 'FontSize',16);
%     legend("\mu = "+ mus(1), "\mu = " + mus(2), ...
%         "Benveniste (\mu_{initial} = "+mu_init+", \rho = "+rho+")", "Ang & Farhang (\mu_{initial} = "+mu_init+", \rho = "+rho+", \alpha = "+alpha_ang+")", ...
%         "Matthews & Xie (\mu_{initial} = "+mu_init+", \rho = "+rho+")", 'FontSize',10, 'Location','southeast');
