function X0 = generate_latent_positions(n, design_type, varargin)
% GENERATE_LATENT_POSITIONS Create latent position matrix with controlled eigenvalues
%
% Inputs:
%   n           - Number of vertices
%   design_type - Type of design:
%                 'optimized': Optimized design (λ3=20, cond#=33) (RECOMMENDED)
%                 'current': Original design (0.3t+0.5, 0.15sin, 0.1cos)
%                 'scaled': Scaled version with larger eigenvalues
%                 'orthogonal': Orthogonalized basis with prescribed eigenvalues
%                 'prescribed': Direct eigenvalue prescription (requires seed)
%   varargin    - Optional parameters:
%                 'seed', value: Random seed for 'prescribed' design
%                 'eigenvalues', [λ1, λ2, λ3]: Target eigenvalues for 'orthogonal'
%
% Output:
%   X0 - n x 3 latent position matrix
%
% Example:
%   X0 = generate_latent_positions(1000, 'optimized');  % Recommended
%   X0 = generate_latent_positions(1000, 'current');
%   X0 = generate_latent_positions(1000, 'orthogonal', 'eigenvalues', [400, 250, 100]);
%   X0 = generate_latent_positions(1000, 'prescribed', 'seed', 42);

% Parse optional arguments
p = inputParser;
addParameter(p, 'seed', 42);
addParameter(p, 'eigenvalues', [400, 250, 100]);
parse(p, varargin{:});

seed = p.Results.seed;
target_eigs = p.Results.eigenvalues;

t = (1:n)' / n;

switch lower(design_type)
    case 'optimized'
        % Optimized design: λ3=20, λ2-λ3=26.3, condition#=33.3
        % Edge probabilities in [0.39, 0.95]
        % Found via optimize_balanced_design.m (Design 9)
        X0 = [0.42*t + 0.46, ...
              0.27 * sin(2*pi*t) + 0.46, ...
              0.20 * cos(4*pi*t)];

    case 'current'
        % Original design from prior simulations
        % λ3=5, condition#=157 (poor numerical properties)
        X0 = [0.3*t + 0.5, ...
              0.15 * sin(2*pi*t) + 0.6, ...
              0.1 * cos(4*pi*t)];

    case 'scaled'
        % Scaled version with larger coefficients for better eigenvalue separation
        % This roughly doubles all scaling factors
        X0 = [0.5*t + 0.6, ...
              0.3 * sin(2*pi*t) + 0.6, ...
              0.2 * cos(4*pi*t)];

    case 'orthogonal'
        % Orthogonalized basis with prescribed eigenvalues
        X_raw = [t, sin(2*pi*t), cos(4*pi*t)];
        [Q, ~] = qr(X_raw, 0);  % Orthonormal columns

        % Scale to get desired eigenvalues
        scaling = sqrt(target_eigs);
        X0 = Q * diag(scaling);

        % Shift to positive range for edge probabilities
        X0 = X0 + [0.6, 0.6, 0];

    case 'prescribed'
        % Direct eigenvalue prescription with random rotation
        rng(seed);
        [U, ~] = qr(randn(n, 3));
        Lambda = diag(target_eigs);
        X0 = U * sqrt(Lambda);

        % Shift to positive range
        X0 = X0 + [0.6, 0.6, 0];

    otherwise
        error('Unknown design_type: %s. Choose: current, scaled, orthogonal, prescribed', design_type);
end

% Validate edge probabilities are reasonable
S = diag([1, 1, -1]);
P = X0 * S * X0';
if min(P(:)) < -0.1 || max(P(:)) > 1.1
    warning('Edge probabilities out of range [0,1]: [%.3f, %.3f]', min(P(:)), max(P(:)));
end

end
