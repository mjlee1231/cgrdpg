function [X_opt, Z_opt, fval, exitflag, output] = fit_grdpg_fminunc_surrogate(A, B, d, p, tau, options)
% FIT_GRDPG_FMINUNC_SURROGATE Fit GRDPG with surrogate likelihood (matching R implementation)
%
% Implements surrogate/majorization algorithm with outer iterations:
%   1. Initialize with ASE and estimate signature matrix S from eigenvalue signs
%   2. Outer loop (max 30 iterations, convergence on max row change < 0.01):
%      a. Fix Y_hat = X * S and Z_hat = B * X * (X'X)^{-1}
%      b. Optimize X with FIXED Y_hat and Z_hat using fminunc
%      c. Check convergence: max row change < tol
%
% Objective (with Y_hat, Z_hat FIXED):
%   ell(X; Y_hat, Z_hat) = sum((A - X*Y_hat') .* psi(X*Y_hat') + Psi(X*Y_hat'))
%                          - 0.5 * ||B - Z_hat*X'||_F^2
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
    % Default options for inner fminunc (trust-region, requires gradient)
    options = optimoptions('fminunc', ...
        'Algorithm', 'trust-region', ...
        'Display', 'off', ...
        'MaxIterations', 100, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-10, ...
        'SpecifyObjectiveGradient', true, ...
        'HessianApproximation', 'lbfgs');
end

n = size(A, 1);
p_cov = size(B, 1);
q = d - p;

% Step 0: Initialize with ASE and estimate signature matrix
fprintf('Initializing with ASE...\n');
[X_current, S_estimated] = initialize_ase(A, d);

fprintf('  Estimated signature: S = diag([');
fprintf('%+d ', diag(S_estimated)');
fprintf('])\n');

% Outer loop parameters
max_outer_iter = 30;
tol_outer = 0.01;  % Convergence tolerance on max row change

% Outer loop: Surrogate optimization
fprintf('\nStarting surrogate optimization with outer iterations...\n');
fprintf('  Parameters: n=%d, p_cov=%d, d=%d, tau=%.6f\n', n, p_cov, d, tau);
fprintf('  Outer iterations: max %d, tol=%.4f (max row change)\n\n', max_outer_iter, tol_outer);

fprintf('%-5s | %-15s | %-15s | %-10s\n', 'Iter', 'Objective', 'Max Row Change', 'Inner Its');
fprintf('%s\n', repmat('-', 1, 60));

for outer_iter = 1:max_outer_iter
    % Fix Y_hat and Z_hat for this iteration
    Y_hat = X_current * S_estimated;
    % Stably solve Z_hat from: Z_hat * X_current' = B
    % Using backslash operator to avoid ill-conditioning from X'X
    Z_hat = (X_current \ B')';

    % Inner optimization: minimize over X with FIXED Y_hat, Z_hat
    x0 = X_current(:);
    objective = @(x) surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau);

    [x_opt, fval, inner_exitflag, inner_output] = fminunc(objective, x0, options);
    X_new = reshape(x_opt, n, d);

    % Check convergence: max row change
    row_changes = sqrt(sum((X_new - X_current).^2, 2));
    max_row_change = max(row_changes);

    fprintf('%5d | %+15.6e | %15.6e | %10d\n', ...
        outer_iter, -fval, max_row_change, inner_output.iterations);

    % Convergence check
    if max_row_change < tol_outer
        fprintf('\nConverged: max row change (%.6e) < tol (%.6e)\n', max_row_change, tol_outer);
        exitflag = 1;
        break;
    end

    % Update for next iteration
    X_current = X_new;
end

% Check if reached max iterations without convergence
if outer_iter == max_outer_iter && max_row_change >= tol_outer
    fprintf('\nReached max outer iterations (%d) without convergence\n', max_outer_iter);
    exitflag = 0;
end

% Final outputs
X_opt = X_current;
% Stably solve Z_opt from: Z_opt * X_opt' = B
Z_opt = (X_opt \ B')';

% Create output structure
output.iterations = outer_iter;
output.funcCount = outer_iter * inner_output.funcCount;
output.firstorderopt = max_row_change;  % Use max row change as optimality measure
output.algorithm = 'Surrogate with fminunc (trust-region)';

fprintf('\nOptimization complete:\n');
fprintf('  Exit flag: %d\n', exitflag);
fprintf('  Outer iterations: %d\n', outer_iter);
fprintf('  Final max row change: %.6e\n', max_row_change);
fprintf('  Final objective: %.6e\n', -fval);

end

function [X_init, S_estimated] = initialize_ase(A, d)
    % Initialize with Adjacency Spectral Embedding and estimate signature
    n = size(A, 1);

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

    % Estimate signature from top d eigenvalues
    S_estimated = diag(sign(eigvals(1:d)));

    % Take top d eigenvectors (unsigned version: U |Lambda|^{1/2})
    X_init = V(:, 1:d) * diag(sqrt(abs(eigvals(1:d))));
end

function [f, g] = surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau)
    % Compute surrogate objective and gradient with FIXED Y_hat and Z_hat
    %
    % Objective (Y_hat and Z_hat are FIXED parameters, NOT functions of X):
    %   f = sum((A - X*Y_hat') .* psi(X*Y_hat') + Psi(X*Y_hat'))
    %       - 0.5 * ||B - Z_hat*X'||_F^2
    %
    % Gradient:
    %   grad_X_net = 2 * W * Y_hat  where W = (A - S_mat) .* dpsi(S_mat)
    %   grad_X_cov = (B - Z_hat*X')' * Z_hat
    %   grad_X = grad_X_net + grad_X_cov (for MAXIMIZING)
    %   For MINIMIZING: negate everything

    % Unpack X
    X = reshape(x, n, d);

    % Network probabilities: S_mat(i,j) = x_i^T * y_hat_j
    S_mat = X * Y_hat';
    % Set diagonal to zero (no self-loops)
    S_mat(1:n+1:end) = 0;

    % Compute psi and Psi values
    [psi_val, Psi_val, dpsi_val] = psi_functions(S_mat, tau);

    % Network component: sum((A - S) .* psi(S) + Psi(S))
    net_obj = sum((A(:) - S_mat(:)) .* psi_val(:) + Psi_val(:));

    % Covariate component: -0.5 * ||B - Z_hat*X'||_F^2
    B_pred = Z_hat * X';
    cov_obj = -0.5 * sum((B(:) - B_pred(:)).^2);

    % Total objective (we MINIMIZE, R code MAXIMIZES, so negate)
    f = -(net_obj + cov_obj);

    % Compute gradient if requested
    if nargout > 1
        % 1. Compute full weight matrix (including diagonal)
        W_full = (A - S_mat) .* dpsi_val;

        % 2. Network gradient excluding diagonal contribution
        % Match R's single node sweep by accurately excluding self-loops
        grad_X_net = -W_full * Y_hat;

        % 3. Correction: Remove diagonal contribution that was incorrectly included
        % After matrix multiplication, diagonal terms (i=j) should not contribute
        diag_W = diag(W_full);
        grad_X_net = grad_X_net + diag_W .* Y_hat;

        % 4. Covariate gradient
        % For MAXIMIZING: grad_X_cov = (B - Z_hat*X')' * Z_hat
        % For MINIMIZING: negate to get -(B - Z_hat*X')' * Z_hat
        resid_cov = B - B_pred;
        grad_X_cov = -resid_cov' * Z_hat;

        % 5. Total gradient for X
        grad_X = grad_X_net + grad_X_cov;

        % Pack gradient
        g = grad_X(:);
    end
end
