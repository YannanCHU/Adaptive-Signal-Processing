% three GASS (gradient adaptive step-size) algorithms

clc;
clear;
close all;

% number of samples
N = 2000;
% step sizes
mus = [0.01, 0.1];
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

% store the weight error - Wo - W(n)
weighErrorLMS = cell(length(mus), numOfRealisations);
averageWeightErrorLMS = cell(length(mus),1);

errorLMS_sq = cell(length(mus), numOfRealisations);
averageErrorLMS = cell(length(mus),1);

for i = 1:length(mus)
    for j = 1:numOfRealisations
        % the input signal is the noise
        % thus, the real order of MA should + 1
        [~, errorLMS, weights] = LMS_arma(noise1(:,j), x100(:,j), mus(i), 0, MA_order+1);
        errorLMS_sq{i,j} = (errorLMS .^ 2);
        weighErrorLMS{i,j} = a1 - weights(2,:);
    end
    averageErrorLMS{i,1} = mean(cat(3, errorLMS_sq{i,:}),3);
    averageWeightErrorLMS{i,1} = mean(cat(3, weighErrorLMS{i,:}),3);
end

%% adaptive step size algorithm
errors_ben = cell(1,numOfRealisations);
errors_ang = cell(1,numOfRealisations);
errors_mat = cell(1,numOfRealisations);
weightsError_ben = cell(1,numOfRealisations);
weightsError_ang = cell(1,numOfRealisations);
weightsError_mat = cell(1,numOfRealisations);

mu_inits = [0.01, 0.3];  % initial guess of the step-size
rho = 0.005;     % learning rate for adaptive step-size algorithm
alpha_ang = 0.7;% the additional parameter of Ang-Farhang algorithm
mus1 = zeros(numOfRealisations, N+1);
mus2 = zeros(numOfRealisations, N+1);
mus3 = zeros(numOfRealisations, N+1);

for k = 1:length(mu_inits)
    mu_init = mu_inits(k);
    for j = 1:numOfRealisations
        [~, errors, weights, mus1(j,:)] = GASS_arma(noise1(:,j), x100(:,j), rho, mu_init, 0, MA_order+1, 'Benveniste', []);
        weightsError_ben{j} = (a1 - weights(2,:)).';
        errors_ben{j} = (errors .^ 2);

        [~, errors, weights, mus2(j,:)] = GASS_arma(noise1(:,j), x100(:,j), rho, mu_init, 0, MA_order+1, 'AngFarhang', alpha_ang);
        weightsError_ang{j} = (a1 - weights(2,:)).';
        errors_ang{j} = (errors .^ 2);

        [~, errors, weights, mus3(j,:)] = GASS_arma(noise1(:,j), x100(:,j), rho, mu_init, 0, MA_order+1, 'MatthewsXie', []);
        weightsError_mat{j} = (a1 - weights(2,:)).';
        errors_mat{j} = (errors .^ 2);
    end

    mus1_avg = mean(mus1,1);
    mus2_avg = mean(mus2,1);
    mus3_avg = mean(mus3,1);

    averageErrors_ben = mean(cell2mat(errors_ben),2);
    averageErrors_ang = mean(cell2mat(errors_ang),2);
    averageErrors_mat = mean(cell2mat(errors_mat),2);
    averageWeightsError_ben = mean(cell2mat(weightsError_ben),2);
    averageWeightsError_ang = mean(cell2mat(weightsError_ang),2);
    averageWeightsError_mat = mean(cell2mat(weightsError_mat),2);

    figure(1);
    subplot(length(mu_inits),1,k);
    plot(1:N, averageWeightErrorLMS{1,1}, 1:N, averageWeightErrorLMS{2,1}, ...
        1:N, averageWeightsError_ben, 1:N, averageWeightsError_ang, ...
        1:N, averageWeightsError_mat, ...
        'LineWidth', 2);
    xlabel("Sample index", 'FontSize',14);
    ylabel("Weight Error", 'FontSize', 14);
    title("The weight error curves for standard LMS and GASS (\mu_{initial} in GASS is "+mu_init+")", 'FontSize',16);
    legend("\mu = "+ mus(1), "\mu = " + mus(2), ...
        "Benveniste (\mu_{initial} = "+mu_init+", \rho = "+rho+")", "Ang & Farhang (\mu_{initial} = "+mu_init+", \rho = "+rho+", \alpha = "+alpha_ang+")", ...
        "Matthews & Xie (\mu_{initial} = "+mu_init+", \rho = "+rho+")", 'FontSize',12);
     xlim([0,250]);

    figure(2);
    subplot(length(mu_inits),1,k);
    plot(1:N, pow2db(abs(averageWeightErrorLMS{1,1})), 1:N, pow2db(abs(averageWeightErrorLMS{2,1})), ...
        1:N, pow2db(abs(averageWeightsError_ben)), 1:N, pow2db(abs(averageWeightsError_ang)), ...
        1:N, pow2db(abs(averageWeightsError_mat)), ...
        'LineWidth', 2);
    xlabel("Sample index", 'FontSize',14);
    ylabel("Weight Error in dB", 'FontSize', 14);
    title("The weight error curves for standard LMS and GASS (\mu_{initial} in GASS is "+mu_init+")", 'FontSize',16);
    legend("\mu = "+ mus(1), "\mu = " + mus(2), ...
        "Benveniste (\mu_{initial} = "+mu_init+", \rho = "+rho+")", "Ang & Farhang (\mu_{initial} = "+mu_init+", \rho = "+rho+", \alpha = "+alpha_ang+")", ...
        "Matthews & Xie (\mu_{initial} = "+mu_init+", \rho = "+rho+")", 'FontSize',12);
    xlim([0,1000]);

    figure(3);
    subplot(length(mu_inits),1,k);
    plot(1:N, pow2db(averageErrorLMS{1,1}), 1:N, pow2db(averageErrorLMS{2,1}), ...
        1:N, pow2db(averageErrors_ben), 1:N, pow2db(averageErrors_ang), ...
        1:N, pow2db(averageErrors_mat), ...
        'LineWidth', 2);
    xlabel("Sample index", 'FontSize',14);
    ylabel("Squared prediction error in dB", 'FontSize', 14);
    title("The squared prediction error curves for standard LMS and GASS (\mu_{initial} in GASS is "+mu_init+")", 'FontSize',16);
    legend("\mu = "+ mus(1), "\mu = " + mus(2), ...
        "Benveniste (\mu_{initial} = "+mu_init+", \rho = "+rho+")", "Ang & Farhang (\mu_{initial} = "+mu_init+", \rho = "+rho+", \alpha = "+alpha_ang+")", ...
        "Matthews & Xie (\mu_{initial} = "+mu_init+", \rho = "+rho+")", 'FontSize',11, 'Location','southeast');
    xlim([0,2000]);
end