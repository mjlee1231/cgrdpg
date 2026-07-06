function [X_opt, Z_opt, fval, exitflag, output] = fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options)
% FIT_GRDPG_FMINUNC_SURROGATE Fit GRDPG with surrogate likelihood (matching R implementation)
%
% Uses the same surrogate objective as the R cgrdpg package:
%   Objective = sum((A - S) .* psi(S) + Psi(S)) - 0.5 * ||B - Z*X^T||^2
% where S = X * sign_diag * X^T and psi/Psi are smoothed log-likelihood functions
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
%   fval - final objective value (surrogate negative log-likelihood)
%   exitflag - optimization exit flag
%   output - optimization output structure

if nargin < 5 || isempty(tau)
    tau = 0.001;
end

if nargin < 6
    % Default options for fminunc
    % TESTING: Using numerical gradients instead of analytical
    options = optimoptions('fminunc', ...
        'Algorithm', 'quasi-newton', ...
        'Display', 'iter', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 10000, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'SpecifyObjectiveGradient', false);  % Let MATLAB compute gradient numerically
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

fprintf('Starting fminunc optimization with surrogate objective...\n');
fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
fprintf('  Total parameters: %d\n', length(x0));

% Define objective function with gradient
objective = @(x) surrogate_objective_gradient(x, A, B, S, n, d, p_cov, tau);

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
    eigvals = diag(D);

    % Sort by MAGNITUDE (absolute value) to match R's eigs_sym(..., which="LM")
    % This is critical for GRDPG with negative eigenvalues!
    [~, idx] = sort(abs(eigvals), 'descend');
    eigvals = eigvals(idx);
    V = V(:, idx);

    % Take top d eigenvectors (unsigned version: U |Lambda|^{1/2})
    X_init = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));

    % Initialize Z via least squares
    Z_init = B * X_init / (X_init' * X_init);
end

function [f, g] = surrogate_objective_gradient(x, A, B, S, n, d, p_cov, tau)
    % Compute surrogate objective and gradient matching R implementation
    %
    % Objective:
    %   f = sum((A - S_mat) .* psi(S_mat) + Psi(S_mat)) - 0.5 * ||B - Z*X^T||_F^2
    % where S_mat(i,j) = x_i^T * sign_diag * x_j

    % Unpack parameters
    X = reshape(x(1:n*d), n, d);
    Z = reshape(x(n*d+1:end), p_cov, d);

    % Compute Y = X * S (signed latent positions)
    Y = X * S;

    % Network probabilities: S_mat(i,j) = x_i^T * sign_diag * x_j
    S_mat = X * (S * X');  % Equivalent to X * S * X'
    % Set diagonal to zero (no self-loops)
    S_mat(1:n+1:end) = 0;

    % Compute psi and Psi values
    [psi_val, Psi_val, dpsi_val] = psi_functions(S_mat, tau);

    % Network component: sum((A - S) .* psi(S) + Psi(S))
    % R sums over ALL pairs (i,j) with i≠j, which counts each edge twice for undirected graphs
    % We must match R's approach exactly
    net_obj = sum((A(:) - S_mat(:)) .* psi_val(:) + Psi_val(:));

    % Covariate component: -0.5 * ||B - Z*X^T||_F^2
    B_pred = Z * X';
    cov_obj = -0.5 * sum((B(:) - B_pred(:)).^2);

    % Total objective (we MINIMIZE, R code MAXIMIZES, so negate)
    f = -(net_obj + cov_obj);

    % COMMENTED OUT: Analytical gradient (testing with numerical gradients)
    % Compute gradient if requested
%     if nargout > 1
%         % Gradient w.r.t. X
%         % The derivative of [(A_ij - S_ij) * psi(S_ij) + Psi(S_ij)] w.r.t. S_ij is:
%         % d/dS_ij = -psi(S_ij) + (A_ij - S_ij)*dpsi(S_ij) + psi(S_ij) = (A_ij - S_ij)*dpsi(S_ij)
%         %
%         % Objective sums over ALL (i,j) pairs with i≠j. Derivative w.r.t. x_i includes:
%         % - Terms where i is first index: sum_j [...] * dS_ij/dx_i
%         % - Terms where i is second index: sum_k [...] * dS_ki/dx_i
%         %
%         % Since S and A are symmetric, factor of 2:
%         % grad_X_i = 2 * sum_{j≠i} (A_ij - S_ij)*dpsi(S_ij) * sign_diag * x_j
%         %
%         % In matrix form with W = (A - S) .* dpsi(S) and diagonal = 0:
%         % grad_X = 2 * W * Y * sign_diag (maximizing)
%         % But we're MINIMIZING -f, so negate: grad_X = -2 * W * Y * S
%
%         % Compute weight matrix
%         W_net = (A - S_mat) .* dpsi_val;  % Note: psi terms cancel in derivative!
%         W_net(1:n+1:end) = 0;  % Zero diagonal (no self-loops)
%
%         % Network gradient (negate for minimization)
%         % Since objective sums over ALL pairs (i,j), gradient already includes both directions
%         % NO factor of 2 needed (would be double-counting for undirected graph)
%         grad_X_net = -W_net * Y;
%
%         % Covariate gradient (negate for minimization)
%         resid_cov = B - B_pred;
%         grad_X_cov = -resid_cov' * Z;
%
%         % Total gradient for X (both components already negated)
%         grad_X = grad_X_net + grad_X_cov;
%
%         % Gradient w.r.t. Z (negate for minimization)
%         % For maximizing: d/dZ[-0.5*||B - Z*X^T||^2] = (B - Z*X') * X
%         % For minimizing: negate to get -(B - Z*X') * X
%         grad_Z = -resid_cov * X;
%
%         % Pack gradient
%         g = [grad_X(:); grad_Z(:)];
%     end
end
