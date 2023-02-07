% this function was initially created for Task-3.1-(a)

function [x_hat, errors, weights_h, weights_g] = ACLMS_arma(x_in, x, mu, AR_order, MA_order)
% apply the LMS adaptive prediction
% inputs:
% x_in (1 x N vector): the input signal with N samples
% x (1 x N vector): the reference signal with N samples
% mu (scalar): the specified step size
% AR_order (scalar):   order of AR model
% MA_order (scalar):   order of MA model
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector

N = length(x);  % signal length
weights_h = zeros(AR_order+MA_order, N+1);   % store the weights computed in each iteration
weights_g = zeros(AR_order+MA_order, N+1);   % store the weights computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
in_vector = zeros(AR_order+MA_order, N);    % store past real signals
x_hat = zeros(N,1);             % store the estimated siganl

for delayIndex = 1:1:AR_order
    % slice the reference signal
    in_vector(delayIndex, 1+delayIndex:end) = x(1,1:N-delayIndex);
end

for delayIndex = 1:1:MA_order
    % slice the reference signal
    in_vector(AR_order+delayIndex, delayIndex:end) = x_in(1,1:N-delayIndex+1);
end

for i = 1:1:N
    % estimate the current signal by current weight vector and past signals
    x_hat(i) = weights_h(:,i)' * in_vector(:,i) + weights_g(:,i)' * conj(in_vector(:,i));
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
    weights_h(:,i+1) = weights_h(:,i) + mu * conj(errors(i)) * in_vector(:,i);
    weights_g(:,i+1) = weights_g(:,i) + mu * conj(errors(i)) * conj(in_vector(:,i));
end

% the last weight vector would be never used
weights_h(:,end) = [];
weights_g(:,end) = [];
end