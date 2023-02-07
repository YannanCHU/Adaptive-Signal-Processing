clc;
clear;
close all;

%% voltage simulation
N = 1000;           % number of samples
fo = 50;            % unit: Hz - system frequency
fs = 1e3;           % unit: Hz - sampling frequency
phi = 0;            % phase shift
delta_b = 0;        % phase distortion
delta_c = 0;        % phase distortion
n = 0:N-1;          % sample index
t = n / fs;         % sampling time
V = 1;              % amplitude, Va = Vb = Vc = V


%% balanced condition
deltas = [0; delta_b; delta_c];
vn_bal = ClarkeVoltageGeneration(V*ones(3,1), fo, fs, n, phi, deltas);
coef_vn_bal = circularityMeasure(vn_bal)

% figure(1);
% scatter(real(vn_bal), imag(vn_bal), '.');
% daspect([1,1,1]);
% title(["Complex balanced voltage", ...
%     "(|\rho| = " + round(coef_vn_bal,3) + ")"],'FontSize',14);
% xlabel("Real Part", 'FontSize',14); 
% ylabel("Imaginary Part", 'FontSize',14);

%% unbalanced condition
delta_b = pi / 3;
delta_c = pi / 4;
deltas = [0; delta_b; delta_c];
Va = V;
Vb = 1.4 * V;
Vc = 0.6 * V;
Vs = [Va; Vb; Vc];
vn_unbal = ClarkeVoltageGeneration(Vs, fo, fs, n, phi, deltas);
coef_vn_unbal = circularityMeasure(vn_unbal);

figure(1);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
scatter(real(vn_bal), imag(vn_bal), 'b', 'filled');
daspect([1,1,1]);
title(["Complex balanced voltage (|\rho| = " + round(coef_vn_bal,3) + ")", ...
    "V_a = 1, V_b = 1, V_c = 1; \Delta_b = 0 rad, \Delta_c = 0 rad"],'FontSize',14);
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);

nexttile;
scatter(real(vn_unbal), imag(vn_unbal), 'r', 'filled');
daspect([1,1,1]);
title(["Complex unbalanced voltage (|\rho| = " + round(coef_vn_unbal,3) + ")", ...
    sprintf("V_a = %.2f, V_b = %.2f, V_c = %.2f; \\Delta_b = %.2f rad, \\Delta_c = %.2f rad", Va, Vb, Vc, delta_b, delta_c)],'FontSize',14);
xlabel("Real Part", 'FontSize',14); 
ylabel("Imaginary Part", 'FontSize',14);

%% estimate the frequency
% balanced voltages
mu = 0.08;
filterOrder = 1;
[~, error_bal_clms, weights_clms] = CLMS_arma([], vn_bal, mu, filterOrder, 0);
fo_bal_clms = abs((fs / (2*pi)) .* atan(-imag(weights_clms) ./ real(weights_clms)));
fprintf("The frequency of balanced voltage estimated by CLMS is %.2f Hz\n", fo_bal_clms(end));

[~, error_bal_aclms, weights_aclms_h, weights_aclms_g] = ACLMS_arma([], vn_bal, mu, filterOrder, 0);
h = weights_aclms_h;
g = weights_aclms_g;
fo_bal_aclms = abs(( fs / (2*pi)) .* atan(sqrt( imag(h).^2 - abs(g).^2 ) ./ real(h) ));
fprintf("The frequency of balanced voltage estimated by ACLMS is %.2f Hz\n", fo_bal_aclms(end));


% unbalanced voltages
[~, error_unbal_clms, weights_clms] = CLMS_arma([], vn_unbal, mu, filterOrder, 0);
fo_unbal_clms = abs((fs / (2*pi)) .* atan(-imag(weights_clms) ./ real(weights_clms)));
fprintf("The frequency of unbalanced voltage estimated by CLMS is %.2f Hz\n", fo_unbal_clms(end));

[~, error_unbal_aclms, weights_aclms_h, weights_aclms_g] = ACLMS_arma([], vn_unbal, mu, filterOrder, 0);
h = weights_aclms_h;
g = weights_aclms_g;
fo_unbal_aclms = abs(( fs / (2*pi)) .* atan(sqrt( imag(h).^2 - abs(g).^2 ) ./ real(h) ));
fprintf("The frequency of balanced voltage estimated by ACLMS is %.2f Hz\n", fo_unbal_aclms(end));

figure(2);
tiledlayout(2,2,'TileSpacing','compact');
nexttile;
plot(n, fo_bal_clms, n, fo_bal_aclms, n, 50*ones(1,N), 'g--', 'LineWidth', 2);
title("Frequency estimates for the balanced system voltages (|\rho| = " + round(coef_vn_bal,3) + ")",'FontSize',14);
xlabel("Sample Index",'FontSize',14);
ylabel("Frequency Estimate",'FontSize',14);
legend("CLMS", "ACLMS", "True frequency", 'FontSize', 12);
ylim([10,90]);

nexttile;
plot(n, fo_unbal_clms, n, fo_unbal_aclms, n, 50*ones(1,N), 'g--', 'LineWidth', 2);
title("Frequency estimates for the unbalanced system voltages (|\rho| = " + round(coef_vn_unbal,3) + ")",'FontSize',14);
xlabel("Sample Index",'FontSize',14);
ylabel("Frequency Estimate",'FontSize',14);
legend("CLMS", "ACLMS", "True frequency", 'FontSize', 12);
ylim([10,90]);


% plot the error
nexttile;
plot(n, pow2db(abs(error_bal_clms).^2), n, pow2db(abs(error_bal_aclms).^2), 'LineWidth', 2);
title("Voltage prediction errors for the balanced system voltages (|\rho| = " + round(coef_vn_bal,3) + ")",'FontSize',14);
xlabel("Sample Index",'FontSize',14);
ylabel("Voltage prediction error in dB",'FontSize',14);
legend("CLMS", "ACLMS", 'FontSize', 12);

% plot the error
nexttile;
plot(n, pow2db(abs(error_unbal_clms).^2), n, pow2db(abs(error_unbal_aclms).^2), 'LineWidth', 2);
title("Voltage prediction errors for the unbalanced system voltages (|\rho| = " + round(coef_vn_unbal,3) + ")",'FontSize',14);
xlabel("Sample Index",'FontSize',14);
ylabel("Voltage prediction error in dB",'FontSize',14);
legend("CLMS", "ACLMS", 'FontSize', 12);
