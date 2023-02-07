clc;
clear;
close all;

%% import the data to the workspace
v_low_wind = load('low-wind.mat');
vL = v_low_wind.v_east + 1i * v_low_wind.v_north;
coef_vL = circularityMeasure(vL);

v_medium_wind = load('medium-wind.mat');
vM = v_medium_wind.v_east + 1i * v_medium_wind.v_north;
coef_vM = circularityMeasure(vM);

v_high_wind = load('high-wind.mat');
vH = v_high_wind.v_east + 1i * v_high_wind.v_north;
coef_vH = circularityMeasure(vH);

%% scatter diagram of the real and imaginary parts of wind signals obtained
% in three regions (low, medium and high dynamics regions)
figure(1);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
scatter(real(vL(:)), imag(vL(:)), 'r.');
daspect([1,1,1]);
title(["Complex wind signal for low dynamics  regime", ...
    "(|\rho| = " + round(coef_vL,3) + ")"],'FontSize',14);
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);

nexttile;
scatter(real(vM(:)), imag(vM(:)), 'g.');
daspect([1,1,1]);
title(["Complex wind signal for medium dynamics regime", ...
    "(|\rho| = " + round(coef_vM,3) + ")"],'FontSize',14);
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);

nexttile;
scatter(real(vH(:)), imag(vH(:)), 'b.');
daspect([1,1,1]);
title(["Complex wind signal for high dynamics regime", ...
    "(|\rho| = " + round(coef_vH,3) + ")"],'FontSize',14);
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);

%% CLMS and ACLMS
mu_low = 0.1;
mu_medium = 0.01;
mu_high = 0.001;
mus = [mu_low, mu_medium, mu_high];

vs = [vL, vM, vH];

filterOrders = 1:16;
errors_CLMS = zeros(3, length(filterOrders));
errors_ACLMS = zeros(3, length(filterOrders));

for j = 1:1:size(vs,2)
    for i = 1:1:length(filterOrders)
        [~, error_clms, ~] = CLMS_arma([], vs(:,j).', mus(j), filterOrders(i), 0);
        errors_CLMS(j,i) = pow2db(mean(abs(error_clms).^2));

        [~, error_aclms, ~, ~] = ACLMS_arma([], vs(:,j).', mus(j), filterOrders(i), 0);
        errors_ACLMS(j,i) = pow2db(mean(abs(error_aclms).^2));
        % figure();
        % plot(1:length(vH), pow2db(abs(error_clms).^2), 1:length(vH), pow2db(abs(error_aclms).^2),'LineWidth', 2);
    end
end

%%
figure(2);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
plot(filterOrders, errors_CLMS(1,:), filterOrders, errors_ACLMS(1,:), 'LineWidth', 2);
xlabel("Filter Order - M", 'FontSize', 14);
ylabel("MSPE in dB", 'FontSize', 14);
title(["MSPE against filter order (Low dynamics)","(Step size \mu = " + mus(1) + ")"], 'FontSize', 14);
legend("CLMS", 'ACLMS', 'FontSize', 12);

nexttile;
plot(filterOrders, errors_CLMS(2,:), filterOrders, errors_ACLMS(2,:), 'LineWidth', 2);
xlabel("Filter Order - M", 'FontSize', 14);
ylabel("MSPE in dB", 'FontSize', 14);
title(["MSPE against filter order (Medium dynamics)","(Step size \mu = " + mus(2) + ")"], 'FontSize', 14);
legend("CLMS", 'ACLMS', 'FontSize', 12, "Location", "southeast");

nexttile;
plot(filterOrders, errors_CLMS(3,:), filterOrders, errors_ACLMS(3,:), 'LineWidth', 2);
xlabel("Filter Order - M", 'FontSize', 14);
ylabel("MSPE in dB", 'FontSize', 14);
title(["MSPE against filter order (High dynamics)","(Step size \mu = " + mus(3) + ")"], 'FontSize', 14);
legend("CLMS", 'ACLMS', 'FontSize', 12);

%% plot the estimated data
N = length(vL);
Ms = [7, 5, 5];
vs_est_clms = zeros(N,3);
vs_est_aclms = zeros(N,3);

for j = 1:1:size(vs,2)
    [vs_est_clms(:,j), ~, ~] = CLMS_arma([], vs(:,j).', mus(j), Ms(j), 0);
    
    [vs_est_aclms(:,j), ~, ~, ~] = ACLMS_arma([], vs(:,j).', mus(j), 4, 0);
end

figure();
tiledlayout(2,3,'TileSpacing','compact');
nexttile;
scatter(real(vL(:)), imag(vL(:)), 'k.'); hold on;
scatter(real(vs_est_clms(:,1)), imag(vs_est_clms(:,1)), 'r.'); hold off;
title("CLMS estimates (Low dynamics) \mu = "+mus(1) +", M = "+ Ms(1) +")");
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("True", "CLMS", 'FontSize',12);

nexttile;
scatter(real(vM(:)), imag(vM(:)), 'k.'); hold on;
scatter(real(vs_est_clms(:,2)), imag(vs_est_clms(:,2)), 'g.'); hold off;
title("CLMS estimates (Medium dynamics) \mu = "+mus(2) +", M = "+ Ms(2) +")");
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("True", "CLMS", 'FontSize',12);

nexttile;
scatter(real(vH(:)), imag(vH(:)), 'k.'); hold on;
scatter(real(vs_est_clms(:,3)), imag(vs_est_clms(:,3)), 'b.'); hold off;
title("CLMS estimates (High dynamics) \mu = "+mus(3) +", M = "+ Ms(3) +")");
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("True", "CLMS", 'FontSize',12);


nexttile;
scatter(real(vL(:)), imag(vL(:)), 'k.'); hold on;
scatter(real(vs_est_aclms(:,1)), imag(vs_est_aclms(:,1)), 'r.'); hold off;
title("ACLMS estimates (Low dynamics) \mu = "+mus(1) +", M = "+ 4 +")");
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("True", "ACLMS", 'FontSize',12);

nexttile;
scatter(real(vM(:)), imag(vM(:)), 'k.'); hold on;
scatter(real(vs_est_aclms(:,2)), imag(vs_est_aclms(:,2)), 'g.'); hold off;
title("ACLMS estimates (Medium dynamics) \mu = "+mus(2) +", M = "+ 4 +")");
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("True", "ACLMS", 'FontSize',12);

nexttile;
scatter(real(vH(:)), imag(vH(:)), 'k.'); hold on;
scatter(real(vs_est_aclms(:,3)), imag(vs_est_aclms(:,3)), 'b.'); hold off;
title("ACLMS estimates (High dynamics) \mu = "+mus(3) +", M = "+ 4 +")");
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);
legend("True", "ACLMS", 'FontSize',12);