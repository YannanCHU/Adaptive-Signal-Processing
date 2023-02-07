% this function was initially created for Task-2.1-(b)

function [x_hat, errors, weights] = LMS_tanh(x_in, x, mu, AR_order, MA_order, a, biasMode, weights_init)
% apply the LMS adaptive prediction
% inputs:
% x_in (N x 1 vector): the input signal with N samples
% x (N x 1 vector): the reference signal with N samples
% mu (scalar): the specified step size
% AR_order (scalar):   order of AR model
% MA_order (scalar):   order of MA model
% a (scalar): scaling factor of the activation function
% biasMode (string):   bias input to our model
% weights_init (AR_order+MA_order x 1 vector): initial guess of the weight vector
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order+MA_order x N+1): the estimated weight vector

x_in = x_in.';  % get the (1 x N) row vector
x = x.';        % get the (1 x N) row vector
N = length(x);  % signal length
weights = zeros(AR_order+MA_order, N+1);   % store the weights computed in each iteration
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

if strcmp(biasMode, "bias")
    in_vector = [ones(1,N); in_vector];
    weights = [zeros(1,N+1); weights];
end

if isempty(weights_init) == false
    weights(:,1) = weights_init;
end

for i = 1:1:N
    % estimate the current signal by current weight vector and past signals
    x_hat(i) = a*tanh(weights(:,i).' * in_vector(:,i));
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
%     weights(:,i+1) = weights(:,i) + mu * errors(i) * in_vector(:,i);
    % if the update is derivated from the new cost function
    weights(:,i+1) = weights(:,i) + mu * errors(i) * a * (1-(tanh(weights(:,i).' * in_vector(:,i)))^2) * in_vector(:,i);
end

% the last weight vector would be never used
% weights(:,end) = [];
end