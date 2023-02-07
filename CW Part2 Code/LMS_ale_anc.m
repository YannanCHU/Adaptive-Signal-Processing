% this function was initially created for Task-2.3-(c)

function [x_hat, errors, weights] = LMS_ale_anc(esp, x, mu, M, delay, mode)
% apply the LMS adaptive prediction for Adaptive Line Enhancer (ALE) and
% Adaptive Noise Canceller (ANC)
% inputs:
% esp (N x 1 vector): the reference noise
% x (N x 1 vector): the noisy signal with N samples
% mu (scalar): the specified step size
% M (scalar):   order of filter
% delay (scalar):   the delay
% mode (string):    specify the mode of this function
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector
esp = esp.';    % get the (1 x N) row vector
x = x.';        % get the (1 x N) row vector
N = length(x);  % signal length
weights = zeros(M, N+1);   % store the weights computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
in_vector = zeros(M, N);    % store past real signals
x_hat = zeros(N,1);             % store the estimated siganl

if strcmp(mode, 'ALE')
    for orderIndex = 1:1:M
        % slice the reference signal
        in_vector(orderIndex, 1+((orderIndex-1)+delay):end) = x(1,1:N-(orderIndex-1)-delay);
    end
elseif strcmp(mode, 'ANC')
    for orderIndex = 1:1:M
        % slice the reference signal
        in_vector(orderIndex, 1+(orderIndex-1):end) = esp(1,1:N-(orderIndex-1));
    end
else
    disp("Error, no matched mode. Please specify 'ALE' or 'ANC'.");
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