% this function was initially created for Task-2.2-(a)

function [x_hat, errors, weights, mus] = GASS_arma(x_in, x, rho, mu_init, AR_order, MA_order, methodName, alpha)
% apply the LMS adaptive prediction
% inputs:
% x_in (N x 1 vector): the input signal with N samples
% x (N x 1 vector): the reference signal with N samples
% rho (scalar): the learning rate of adaptive step size
% mu_init (scalar): the initial guess of step size
% AR_order (scalar):   order of AR model
% MA_order (scalar):   order of MA model
% methodName (string): indicate which algorithm is selected
% alpha (scalar): ranges from 0 to 1, required by Ang Rarh
% outputs:
% x_hat (N x 1 vector): the estimated signal
% errors (N x 1 vector): the error computed at each time step
% weights (AR_order x N): the estimated weight vector
% mus (1, N+1): the adaptive step size vector

x_in = x_in.';  % get the (1 x N) row vector
x = x.';        % get the (1 x N) row vector
N = length(x);  % signal length
weights = zeros(AR_order+MA_order, N+1);   % store the weights computed in each iteration
mus = zeros(1, N+1);   % store the step size used in each iteration
mus(1) = mu_init;
Psis = zeros(AR_order+MA_order, N+1);   % store the term ψ computed in each iteration
errors = zeros(N,1);            % store the error computed in each iteration
x_vector = zeros(AR_order+MA_order, N);    % store the input vector
x_hat = zeros(N,1);             % store the estimated siganl

for delayIndex = 1:1:AR_order
    % slice the reference signal
    x_vector(delayIndex, 1+delayIndex:end) = x(1,1:N-delayIndex);
end

for delayIndex = 1:1:MA_order
    % slice the reference signal
    x_vector(AR_order+delayIndex, delayIndex:end) = x_in(1,1:N-delayIndex+1);
end

for i = 1:1:N
    % estimate the current signal by current weight vector and past signals
    x_hat(i) = weights(:,i).' * x_vector(:,i);
    % estimate error between estimated signal and actual signal
    errors(i) = x(i) - x_hat(i);
    % update the weight vector
    weights(:,i+1) = weights(:,i) + mus(i) * errors(i) * x_vector(:,i);

    % update the μ for next iteration
    mus(i+1) = mus(i) + rho*errors(i)*x_vector(:,i).'*Psis(:,i);

    % update the ψ for next iteration
    if strcmp(methodName, 'Benveniste')
        Psis(:,i+1) = ( eye(AR_order+MA_order) - mus(i)*x_vector(:,i)*x_vector(:,i).' ) * Psis(:,i) + errors(i)*x_vector(:,i);
    elseif strcmp(methodName, 'AngFarhang')
        Psis(:,i+1) = alpha * Psis(:,i) + errors(i) * x_vector(:,i);
    elseif strcmp(methodName, 'MatthewsXie')
        Psis(:,i+1) = errors(i) * x_vector(:,i);
    end
end

% the last weight vector would be never used
weights(:,end) = [];
end