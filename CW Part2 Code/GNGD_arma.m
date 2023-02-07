% this function was initially created for Task-2.2-(c)

function [x_hat, errors, weights] = GNGD_arma(x_in, x, mu, rho, AR_order, MA_order)
% apply the generalized normalised gradient descent (GNGD) algorithm
% inputs:
% x_in (N x 1 vector): the input signal with N samples
% x (N x 1 vector): the reference signal with N samples
% mu (scalar): the specified step size
% rho (scalar): the paramater used in regularization factor update
% AR_order (scalar):   order of AR model
% MA_order (scalar):   order of MA model
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector

x_in = x_in.';  % get the (1 x N) row vector
x = x.';        % get the (1 x N) row vector
N = length(x);  % signal length
weights = zeros(AR_order+MA_order, N+1);   % store the weights computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
in_vector = zeros(AR_order+MA_order, N);    % store past real signals
x_hat = zeros(N,1);             % store the estimated siganl
epsilon = ones(N+1,1) / mu;     % set the initial two regularization as 1 / mu
beta = mu;

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
    x_hat(i) = weights(:,i).' * in_vector(:,i);
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
    weights(:,i+1) = weights(:,i) + beta * errors(i) * in_vector(:,i) / (epsilon(i) + in_vector(:,i).'*in_vector(:,i));
    % update the regularization factor
    if i == 1
        continue;
    end
    epsilon(i+1) = epsilon(i) - rho*mu* (errors(i)*errors(i-1)*in_vector(:,i).'*in_vector(:,i-1)) / (epsilon(i-1) + in_vector(:,i-1).'*in_vector(:,i-1))^2;
end

% the last weight vector would be never used
weights(:,end) = [];
end