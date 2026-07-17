function Psi_val = Psi(x, tau)
% PSI Integral of smoothed log-likelihood function psi
%
% Psi(x) = integral of psi(x) = Psi1(x) + Psi2(x)
% Used in surrogate objective computation
%
% Inputs:
%   x   - numeric vector of values
%   tau - smoothing threshold (default: 0.001)
%
% Output:
%   Psi_val - integral values

if nargin < 2
    tau = 0.001;
end

[~, Psi_val, ~, ~] = psi_functions(x, tau);

end
