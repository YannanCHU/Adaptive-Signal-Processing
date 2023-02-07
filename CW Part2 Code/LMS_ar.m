% this function was initially created for Task-2.1-(b)

function [x_hat, errors, weights] = LMS_ar(x, mu, AR_order)
% apply the LMS adaptive prediction
% inputs:
% x (N x 1 vector): the reference signal with N samples
% mu (scalar): the specified step size
% AR_order (scalar);   % order of AR model
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector

x = x.';        % get the (1 x N) row vector
N = length(x);  % signal length
weights = zeros(AR_order, N+1);   % store the weights computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
x_past = zeros(AR_order, N);    % store past real signals
x_hat = zeros(N,1);             % store the estimated siganl

for delayIndex = 1:1:AR_order
    % slice the reference signal
    x_past(delayIndex, 1+delayIndex:end) = x(1,1:N-delayIndex);
end

for i = 1:1:N
    % estimate the current signal by current weight vector and past signals
    x_hat(i) = weights(:,i).' * x_past(:,i);
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
    weights(:,i+1) = weights(:,i) + mu * errors(i) * x_past(:,i);
end

% the last weight vector would be never used
weights(:,end) = [];
end