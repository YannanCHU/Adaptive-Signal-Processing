% this function was initially created for Task-3.3-(c)

function [x_hat, errors, weights] = DFT_CLMS(x, mu, gamma)
% apply the DFT-CLMS algorithm
% inputs:
% x (1 x N vector): the reference signal with N samples
% mu (scalar): the specified step size
% gamma (scalar): the leakage coefficient (ranges from 0 to 1)
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector

N = length(x);  % signal length
weights = zeros(N, N+1);   % store the weights computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
x_hat = zeros(N,1);             % store the estimated siganl

DFT_mat = (1/N) * exp(1i * 2*pi * (0:N-1).' * (0:N-1) / N);

for i = 1:1:N
    % estimate the current signal by current weight vector and past signals
    x_hat(i) = weights(:,i)' * DFT_mat(:,i);
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
    weights(:,i+1) = (1-mu*gamma) * weights(:,i) + mu * conj(errors(i)) * DFT_mat(:,i);
end

% the last weight vector would be never used
weights(:,end) = [];
end