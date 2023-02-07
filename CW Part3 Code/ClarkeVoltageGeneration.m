function vn = ClarkeVoltageGeneration(Vs, fo, fs, n, phi, deltas)
% inputs
% Vs (3 x 1 vector): peak voltages of three sinusoids
% fo (scalar): unit: Hz - system frequency
% fs (scalar): unit: Hz - sampling frequency
% n (1 x N vector):  sample index
% phi (scalar): phase shift
% deltas (3 x 1 vector): phase distortions (the first value = 0)
% output
% vn (1 x N vector): Complex Clark Voltage
delta_b = deltas(2);
delta_c = deltas(3);
va = Vs(1) * cos(2*pi*fo*n / fs + phi);
vb = Vs(2) * cos(2*pi*fo*n / fs + phi + delta_b - 2*pi/3);
vc = Vs(3) * cos(2*pi*fo*n / fs + phi + delta_c + 2*pi/3);

% Clarke Matrix
C = sqrt(2/3) * [sqrt(1/2), sqrt(1/2), sqrt(1/2); ...
                 1,         -0.5,      -0.5; ...
                 0,         sqrt(3)/2, -sqrt(3)/2];

Vtemp = C * [va; vb; vc];
V0 = Vtemp(1,:);
Valpha = Vtemp(2,:);
Vbeta = Vtemp(3,:);

% a complex Clarke voltage
vn = Valpha + 1i * Vbeta;
end
