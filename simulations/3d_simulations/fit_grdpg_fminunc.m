function [X_opt, Z_opt, fval, exitflag, output] = fit_grdpg_fminunc(A, B, d, p, tau, options)
% FIT_GRDPG_FMINUNC Fit GRDPG with covariates using MATLAB's fminunc
%
% Inputs:
%   A - (n x n) adjacency matrix
%   B - (p_cov x n) covariate matrix (features x nodes)
%   d - embedding dimension
%   p - number of positive signature directions (q = d - p negatives)
%   tau - smoothing parameter (default: 0.001)
%   options - fminunc options structure (optional)
%
% Outputs:
%   X_opt - (n x d) optimized latent positions
%   Z_opt - (p_cov x d) optimized covariate coefficients
%   fval - final objective value (negative log-likelihood)
%   exitflag - optimization exit flag
%   output - optimization output structure
%
% Uses MATLAB's fminunc with gradient-based optimization

if nargin < 5 || isempty(tau)
    tau = 0.001;
end

if nargin < 6
    % Default options for fminunc
    options = optimoptions('fminunc', ...
        'Algorithm', 'quasi-newton', ...  % or 'trust-region' if gradient provided
        'Display', 'iter', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 10000, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'SpecifyObjectiveGradient', true);  % We'll provide gradients
end

n = size(A, 1);
p_cov = size(B, 1);
q = d - p;

% Initialize with ASE
fprintf('Initializing with ASE...\n');
[X_init, Z_init] = initialize_ase(A, B, d, p);

% Signature matrix
S = diag([ones(p, 1); -ones(q, 1)]);

% Pack parameters into a vector for fminunc
% Parameters: [X(:); Z(:)]
x0 = [X_init(:); Z_init(:)];

fprintf('Starting fminunc optimization...\n');
fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
fprintf('  Total parameters: %d\n', length(x0));

% Define objective function with gradient
objective = @(x) grdpg_objective_gradient(x, A, B, S, n, d, p_cov, tau);

% Run fminunc
[x_opt, fval, exitflag, output] = fminunc(objective, x0, options);

% Unpack optimized parameters
X_opt = reshape(x_opt(1:n*d), n, d);
Z_opt = reshape(x_opt(n*d+1:end), p_cov, d);

fprintf('\nOptimization complete:\n');
fprintf('  Exit flag: %d\n', exitflag);
fprintf('  Iterations: %d\n', output.iterations);
fprintf('  Function evaluations: %d\n', output.funcCount);
fprintf('  Final objective: %.6f\n', fval);
fprintf('  First-order optimality: %.6e\n', output.firstorderopt);

end

function [X_init, Z_init] = initialize_ase(A, B, d, p)
    % Initialize with Adjacency Spectral Embedding
    n = size(A, 1);
    p_cov = size(B, 1);

    % Augmented adjacency for better initialization
    A_aug = A;
    deg = sum(A, 2);
    A_aug(1:n+1:end) = deg / (n - 1);  % Set diagonal

    % Eigendecomposition
    [V, D] = eig(A_aug);
    [eigvals, idx] = sort(diag(D), 'descend');
    V = V(:, idx);

    % Take top d eigenvectors
    X_init = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

    % Initialize Z via least squares
    Z_init = B * X_init / (X_init' * X_init);
end

function [f, g] = grdpg_objective_gradient(x, A, B, S, n, d, p_cov, tau)
    % Compute negative log-likelihood and gradient for GRDPG with covariates

    % Unpack parameters
    X = reshape(x(1:n*d), n, d);
    Z = reshape(x(n*d+1:end), p_cov, d);

    % Compute Y = X * S (signed latent positions)
    Y = X * S;

    % Network component: edge probabilities s_ij = y_i' * y_j
    P_net = Y * Y';  % n x n matrix of probabilities

    % Covariate component: predicted B = Z * X'
    B_pred = Z * X';  % p_cov x n

    % Clamp probabilities to [tau, 1-tau] for numerical stability
    P_clamped = max(min(P_net, 1 - tau), tau);

    % Negative log-likelihood
    % Network term: -sum_{i<j} [A_ij * log(s_ij) + (1-A_ij) * log(1-s_ij)]
    % Only count upper triangle (undirected)
    triu_idx = triu(true(n), 1);  % Upper triangle indices, excluding diagonal

    A_vec = A(triu_idx);
    P_vec = P_clamped(triu_idx);

    nll_network = -sum(A_vec .* log(P_vec) + (1 - A_vec) .* log(1 - P_vec));

    % Covariate term: sum of squared errors
    sse_cov = sum((B(:) - B_pred(:)).^2);

    % Total objective (negative surrogate log-likelihood)
    f = nll_network + 0.5 * sse_cov;

    % Compute gradient if requested
    if nargout > 1
        % Gradient w.r.t. X
        grad_X = zeros(n, d);

        % Network gradient
        % For each i: sum_j [residual_ij * w_ij * y_j]
        % where residual_ij = A_ij - s_ij, w_ij = 1/(s_ij*(1-s_ij))

        % Set diagonal to zero (no self-loops)
        P_clamped(1:n+1:end) = 0.5;  % Avoid division by zero

        residual = A - P_clamped;
        weight = 1 ./ (P_clamped .* (1 - P_clamped));
        weight(1:n+1:end) = 0;  % No self-loops

        % Weighted residuals
        W = residual .* weight;

        % Gradient from network: -Y^T * W for X, then multiply by S
        grad_X_net = -2 * W * Y * S;  % n x d

        % Covariate gradient
        % SSE gradient: -2 * Z^T * (B - Z*X')
        resid_cov = B - B_pred;  % p_cov x n
        grad_X_cov = -2 * resid_cov' * Z;  % n x d

        grad_X = grad_X_net + grad_X_cov;

        % Gradient w.r.t. Z
        % From covariate term only: -2 * (B - Z*X') * X
        grad_Z = -2 * resid_cov * X;  % p_cov x d

        % Pack gradient
        g = [grad_X(:); grad_Z(:)];
    end
end
