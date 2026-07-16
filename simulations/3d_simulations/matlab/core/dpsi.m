function dpsi_val = dpsi(x, tau)
% DPSI First derivative of smoothed log-likelihood function psi
%
% dpsi(x) ≈ 1/x + 1/(1-x) (first derivative of Bernoulli log-likelihood)
% Smoothed with quadratic patches outside [tau, 1-tau]
%
% Inputs:
%   x   - numeric vector of values
%   tau - smoothing threshold (default: 0.001)
%
% Output:
%   dpsi_val - first derivative values

if nargin < 2
    tau = 0.001;
end

dpsi_val = dpsi1(x, tau) - dpsi2(x, tau);

end

function out = dpsi1(x, tau)
    % Derivative of psi1
    a = -1/(2*tau^2);
    b = 2/tau;

    out = zeros(size(x));
    idx = (x < tau);
    out(idx) = 2*a*x(idx) + b;
    out(~idx) = 1 ./ max(x(~idx), eps);
end

function out = dpsi2(x, tau)
    % Derivative of psi2
    s = 1 - tau;
    a = -1/(2*tau^2);
    b = (1 - 2*tau)/tau^2;

    out = zeros(size(x));
    idx = (x > s);
    out(idx) = 2*a*x(idx) + b;
    out(~idx) = -1 ./ max(1 - x(~idx), eps);
end
