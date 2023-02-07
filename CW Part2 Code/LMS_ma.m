% this function was initially created for Task-2.2-(a)

function [x_hat, errors, weights] = LMS_ma(x_in, x, mu, MA_order)
% apply the LMS adaptive prediction
% inputs:
% x (N x 1 vector): the reference signal with N samples
% mu (scalar): the specified step size
% MA_order (scalar);   % order of MA model
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector

x_in = x_in.';  % get the (1 x N) row vector
x = x.';        % get the (1 x N) row vector
N = length(x);  % signal length
weights = zeros(MA_order, N+1);   % store the weights computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
in_vector = zeros(MA_order, N);    % store past real signals
x_hat = zeros(N,1);             % store the estimated siganl

for delayIndex = 1:1:MA_order
    % slice the reference signal
    in_vector(delayIndex, delayIndex:end) = x_in(1,1:N-delayIndex+1);
end

for i = 1:1:N
    % estimate the current signal by current weight vector and past signals
    x_hat(i) = weights(:,i).' * in_vector(:,i);
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
    weights(:,i+1) = weights(:,i) + mu * errors(i) * in_vector(:,i);
end

% the last weight vector would be never used
weights(:,end) = [];
end