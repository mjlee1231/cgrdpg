function X_new = fisher_sweep_X(A, X, Z, B, S, tau, ls_beta, ls_c, ls_max)
% FISHER_SWEEP_X Fisher scoring with coordinate descent for GRDPG
%
% Matches R's fisher_sweep_X function exactly:
% - Coordinate descent: update one vertex at a time
% - Newton direction: solve G * p = S where G is Fisher information
% - Backtracking line search with Armijo condition
% - Y is updated IMMEDIATELY after each vertex (affects next vertices)
%
% Inputs:
%   A       - (n x n) adjacency matrix
%   X       - (n x d) current latent positions
%   Z       - (p_cov x d) covariate coefficient matrix
%   B       - (p_cov x n) covariate matrix
%   S       - (d x d) signature matrix
%   tau     - smoothing parameter (default: 0.001)
%   ls_beta - line search backtracking factor (default: 0.35)
%   ls_c    - Armijo constant (default: 1e-4)
%   ls_max  - max line search iterations (default: 30)
%
% Output:
%   X_new - (n x d) updated latent positions after one sweep

if nargin < 6 || isempty(tau), tau = 0.001; end
if nargin < 7 || isempty(ls_beta), ls_beta = 0.35; end
if nargin < 8 || isempty(ls_c), ls_c = 1e-4; end
if nargin < 9 || isempty(ls_max), ls_max = 30; end

n = size(X, 1);
d = size(X, 2);
p_cov = size(Z, 1);

% Working state that we keep updating as we sweep
X_cur = X;
Y_cur = X_cur * S;
Y_cur_t = Y_cur';
ZtZ = Z' * Z;

% Sweep through all vertices
for i = 1:n
    x_i = X_cur(i, :)';  % Column vector

    % ===== Gradient computation at CURRENT state =====
    % Edge probabilities: s_j = x_i^T * y_j
    s = Y_cur * x_i;  % (n x 1)

    % Weights: w = dpsi(s, tau)
    w = dpsi(s, tau);
    w(i) = 0;  % Exclude self-loop

    % Residuals: r = (A_ij - s_j) * w_j
    r = (A(i, :)' - s) .* w;  % (n x 1)

    % Network score: sum_j r_j * y_j
    s_net = Y_cur_t * r;  % (d x 1)

    % Covariate score: Z^T * (b_i - Z * x_i)
    b_i = B(:, i);
    resid = b_i - Z * x_i;
    s_cov = Z' * resid;  % (d x 1)

    % Total score
    S_score = s_net + s_cov;  % (d x 1)

    % ===== Fisher information =====
    % Clip s for numerical stability in Fisher info only
    s_clipped = max(min(s, 1 - tau), tau);
    w_fisher = dpsi(s_clipped, tau);
    w_fisher(i) = 0;  % Exclude self-loop

    % Network Fisher info: Y^T * diag(w_fisher) * Y
    G_net = Y_cur_t * (Y_cur .* w_fisher);  % (d x d)

    % Total Fisher info
    G = G_net + ZtZ;  % (d x d)

    % ===== Solve G * p = S for Newton direction =====
    % Try Cholesky first (fast, stable if G is SPD)
    [R, flag] = chol(G);
    if flag == 0
        % Cholesky succeeded
        p = R \ (R' \ S_score);
    else
        % Cholesky failed, use general solver
        p = G \ S_score;
    end

    % Check if direction is ascent
    sp = S_score' * p;
    if ~isfinite(sp) || sp <= 0
        % Fall back to gradient
        p = S_score;
        sp = S_score' * p;
        if ~isfinite(sp) || sp <= 0
            % Skip this vertex
            continue;
        end
    end

    % ===== Backtracking line search (Armijo for MAXIMIZATION) =====
    % Condition: f(x_i + eta*p) >= f(x_i) + ls_c * eta * (S^T * p)
    best_eta = 0;
    eta = 1.0;

    % Current objective value
    f0 = surrogate_objective_single(A, X_cur, Z, B, S, tau);

    for bt = 1:ls_max
        % Try step
        x_i_try = x_i + eta * p;
        X_try = X_cur;
        X_try(i, :) = x_i_try';

        % Evaluate objective
        f_try = surrogate_objective_single(A, X_try, Z, B, S, tau);

        % Check Armijo condition
        if f_try >= f0 + ls_c * eta * sp
            best_eta = eta;
            break;
        end

        % Reduce step size
        eta = eta * ls_beta;
    end

    % ===== Accept step and UPDATE working state =====
    X_cur(i, :) = (x_i + best_eta * p)';
    Y_cur(i, :) = X_cur(i, :) * S;  % Update Y IMMEDIATELY
    Y_cur_t = Y_cur';  % Update transpose
end

X_new = X_cur;

end

function f = surrogate_objective_single(A, X, Z, B, S, tau)
% Compute surrogate objective value
% f = sum((A - X*Y') .* psi(X*Y') + Psi(X*Y')) - 0.5 * ||B - Z*X'||_F^2
% where Y = X * S

n = size(X, 1);
Y = X * S;

% Network part
XY = X * Y';  % Edge probabilities (n x n)

% Exclude diagonal (no self-loops)
mask = ~eye(n);
XY_vec = XY(mask);
A_vec = A(mask);

% Compute psi and Psi
[psi_vals, Psi_vals, ~, ~] = psi_functions(XY_vec, tau);

net_part = sum((A_vec - XY_vec) .* psi_vals + Psi_vals);

% Covariate part
resid = B - Z * X';
cov_part = -0.5 * sum(resid(:).^2);

f = net_part + cov_part;

end
