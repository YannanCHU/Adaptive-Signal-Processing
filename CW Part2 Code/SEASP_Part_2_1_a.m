clc;
clear;
close all;

a1 = 0.1;
a2 = 0.8;
noiseVar = 0.25;

A = [1, -a1, -a2; a1, a2-1, 0; a2, a1, -1];
b = [noiseVar; 0; 0];
rs = A \ b;

R = [rs(1), rs(2); rs(2), rs(1)];
[eigVec, eigVal] = eig(R);

lambda_max = 2 / max(diag(eigVal));

disp("The auto-correlation matrix is ");
disp(R);

fprintf("The range of step size should be (0, %.2f)\n", lambda_max);