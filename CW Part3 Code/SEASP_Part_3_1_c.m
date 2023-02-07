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
coef_vn_bal = circularityMeasure(vn_bal);

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
