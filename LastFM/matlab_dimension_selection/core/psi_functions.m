function [psi_val, Psi_val, dpsi_val, ddpsi_val] = psi_functions(x, tau)
% PSI_FUNCTIONS Compute smoothed log-likelihood functions for GRDPG
%
% Inputs:
%   x - probability values (can be outside [0,1])
%   tau - smoothing threshold (default: 0.001)
%
% Outputs:
%   psi_val - psi(x) = psi1(x) - psi2(x) ≈ log(x/(1-x))
%   Psi_val - Psi(x) = integral of psi(x)
%   dpsi_val - derivative of psi(x) ≈ 1/x + 1/(1-x)
%   ddpsi_val - second derivative of psi(x) ≈ -1/x² + 1/(1-x)²
%
% These functions provide a smooth, twice-differentiable extension
% of the Bernoulli log-likelihood outside [0,1] using quadratic patches

if nargin < 2
    tau = 0.001;
end

% psi1(x): smoothed log(x)
% Quadratic extension for x < tau, log(x) for x >= tau
psi1_val = psi1(x, tau);
psi2_val = psi2(x, tau);
psi_val = psi1_val - psi2_val;

% Psi1(x) and Psi2(x): primitives (integrals)
if nargout > 1
    Psi1_val = Psi1(x, tau);
    Psi2_val = Psi2(x, tau);
    Psi_val = Psi1_val - Psi2_val;
end

% dpsi(x): first derivative
if nargout > 2
    dpsi1_val = dpsi1(x, tau);
    dpsi2_val = dpsi2(x, tau);
    dpsi_val = dpsi1_val - dpsi2_val;
end

% ddpsi(x): second derivative
if nargout > 3
    ddpsi1_val = ddpsi1(x, tau);
    ddpsi2_val = ddpsi2(x, tau);
    ddpsi_val = ddpsi1_val - ddpsi2_val;
end

end

%% Helper functions

function out = psi1(x, tau)
    % psi1(x): smoothed log(x)
    % Quadratic for x < tau, log(x) for x >= tau
    a = -1/(2*tau^2);
    b =  2/tau;
    c =  log(tau) - 3/2;

    out = zeros(size(x));
    idx = (x < tau);
    out(idx) = a*x(idx).^2 + b*x(idx) + c;
    out(~idx) = log(max(x(~idx), eps));
end

function out = psi2(x, tau)
    % psi2(x): smoothed log(1-x)
    % log(1-x) for x <= 1-tau, quadratic for x > 1-tau
    s = 1 - tau;
    a = -1/(2*tau^2);
    b = (1 - 2*tau)/tau^2;
    c = log(tau) + ((3*tau - 1)*(1 - tau))/(2*tau^2);

    out = zeros(size(x));
    idx = (x > s);
    out(idx) = a*x(idx).^2 + b*x(idx) + c;
    out(~idx) = log(max(1 - x(~idx), eps));
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

function out = Psi1(x, tau)
    % Primitive (integral) of psi1
    a = -1/(2*tau^2);
    b = 2/tau;
    c = log(tau) - 3/2;

    % Quadratic primitive
    Q = @(t) a/3*t.^3 + b/2*t.^2 + c*t;

    out = zeros(size(x));
    Ctau = Q(tau);  % Continuity constant

    idx = (x < tau);
    out(idx) = Q(x(idx));

    % Integral of log(t) from tau to x: x*log(x) - x - (tau*log(tau) - tau)
    out(~idx) = Ctau + (x(~idx).*log(max(x(~idx), eps)) - x(~idx)) - ...
                (tau*log(tau) - tau);
end

function out = Psi2(x, tau)
    % Primitive (integral) of psi2
    s = 1 - tau;
    a = -1/(2*tau^2);
    b = (1 - 2*tau)/tau^2;
    c = log(tau) + ((3*tau - 1)*(1 - tau))/(2*tau^2);

    % Quadratic primitive
    Q = @(t) a/3*t.^3 + b/2*t.^2 + c*t;

    out = zeros(size(x));
    idx = (x <= s);

    % Integral of log(1-t) from 0 to x: -(1-x)*log(1-x) - x
    out(idx) = -(1 - x(idx)) .* log(max(1 - x(idx), eps)) - x(idx);

    if any(~idx)
        % Constant: integral from 0 to s
        Clog = -(1 - s) * log(tau) - s;  % since 1-s = tau
        out(~idx) = Clog + (Q(x(~idx)) - Q(s));
    end
end

function out = ddpsi1(x, tau)
    % Second derivative of psi1
    % For x < tau: psi1 = a*x^2 + b*x + c, so ddpsi1 = 2*a = -1/tau^2
    % For x >= tau: psi1 = log(x), so ddpsi1 = -1/x^2
    a = -1/(2*tau^2);

    out = zeros(size(x));
    idx = (x < tau);
    out(idx) = 2*a;  % Constant: -1/tau^2
    out(~idx) = -1 ./ max(x(~idx), eps).^2;
end

function out = ddpsi2(x, tau)
    % Second derivative of psi2
    % For x > 1-tau: psi2 = a*x^2 + b*x + c, so ddpsi2 = 2*a = -1/tau^2
    % For x <= 1-tau: psi2 = log(1-x), so ddpsi2 = -1/(1-x)^2
    s = 1 - tau;
    a = -1/(2*tau^2);

    out = zeros(size(x));
    idx = (x > s);
    out(idx) = 2*a;  % Constant: -1/tau^2
    out(~idx) = -1 ./ max(1 - x(~idx), eps).^2;
end
