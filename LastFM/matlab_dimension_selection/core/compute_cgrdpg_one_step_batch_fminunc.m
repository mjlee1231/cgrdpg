function [X_new, Z_new] = compute_cgrdpg_one_step_batch_fminunc(A, B, X_init, S_estimated, tau, options)
% COMPUTE_CGRDPG_ONE_STEP_BATCH_FMINUNC One-step batch surrogate with fminunc
%
% Performs ONE outer iteration of batch surrogate optimization:
%   1. Fix Y_hat = X_init * S and Z_hat = B * X_init * (X_init'X_init)^{-1}
%   2. Optimize X using fminunc with FIXED Y_hat and Z_hat
%   3. Return optimized X (one outer iteration)
%
% This is the TRUE batch method using simultaneous update of all vertices.
%
% Inputs:
%   A            - (n x n) adjacency matrix
%   B            - (p_cov x n) covariate matrix
%   X_init       - (n x d) initial positions (from ASE unsigned)
%   S_estimated  - (d x d) estimated signature matrix
%   tau          - clipping parameter for edge probabilities
%   options      - fminunc options structure (optional)
%
% Outputs:
%   X_new  - (n x d) updated latent positions (one outer iteration)
%   Z_new  - (p_cov x d) updated covariate loadings

if nargin < 6 || isempty(options)
    % Default options for fminunc (trust-region, requires gradient and Hessian)
    options = optimoptions('fminunc', ...
        'Algorithm', 'trust-region', ...
        'Display', 'off', ...
        'MaxIterations', 100, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-10, ...
        'SpecifyObjectiveGradient', true, ...
        'HessianFcn', 'objective');
end

n = size(A, 1);
p_cov = size(B, 1);
d = size(X_init, 2);

% Fix Y_hat and Z_hat for this outer iteration
Y_hat = X_init * S_estimated;
Z_hat = (X_init \ B')';  % Solve Z_hat * X_init' = B

% Inner optimization: minimize over X with FIXED Y_hat, Z_hat
x0 = X_init(:);
objective = @(x) surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau);

% Run fminunc (this does many inner iterations)
[x_opt, ~, ~, ~] = fminunc(objective, x0, options);
X_new = reshape(x_opt, n, d);

% Update Z based on new X
Z_new = (X_new \ B')';

end

%% Surrogate objective function with gradient and Hessian

function [f, g, H] = surrogate_objective_gradient(x, A, B, Y_hat, Z_hat, n, d, tau)
    % Compute surrogate objective, gradient, and Hessian with FIXED Y_hat and Z_hat

    % 1. Unpack X (column-major)
    X = reshape(x, n, d);

    % 2. Network probabilities
    S_mat = X * Y_hat';

    % Create off-diagonal mask (diagonal = 0, off-diagonal = 1)
    is_off_diag = ~eye(n);

    % 3. Compute psi and Psi values
    [psi_val, Psi_val, dpsi_val, ddpsi_val] = psi_functions(S_mat, tau);

    % 4. Network component calculation (exclude diagonal)
    net_obj = sum((A(is_off_diag) - S_mat(is_off_diag)) .* psi_val(is_off_diag) + Psi_val(is_off_diag));

    % 5. Covariate component: -0.5 * ||B - Z_hat*X'||_F^2
    B_pred = Z_hat * X';  % (p_cov x n)
    cov_obj = -0.5 * sum((B - B_pred).^2, 'all');

    % 6. Total objective (fminunc minimization, so negate)
    f = -(net_obj + cov_obj);

    % 7. Compute gradient if requested
    if nargout > 1
        % Network gradient weight matrix
        W_core = (A - S_mat) .* dpsi_val;
        W_core(~is_off_diag) = 0;  % Exclude diagonal

        % Network gradient (MINIMIZING - negate)
        grad_X_net = -W_core * Y_hat;

        % Covariate gradient (MINIMIZING - negate)
        grad_X_cov = -(B - B_pred)' * Z_hat;

        % Total gradient
        grad_X = grad_X_net + grad_X_cov;

        % Pack gradient (column-major)
        g = grad_X(:);
    end

    % 8. Compute Hessian if requested
    if nargout > 2
        % Block diagonal Hessian approximation
        H = zeros(n*d, n*d);

        % Covariate Hessian contribution (constant for all blocks)
        ZtZ = Z_hat' * Z_hat;

        % Network Hessian: different for each vertex
        for i = 1:n
            % Get ddpsi weights for vertex i (exclude diagonal)
            w_i = ddpsi_val(i, :)';
            w_i(i) = 0;  % Exclude self-loop

            % Block Hessian for vertex i
            H_i_net = -Y_hat' * (Y_hat .* w_i);  % (d x d), negative for minimization
            H_i = H_i_net + ZtZ;

            % Place in block diagonal
            idx = (i-1)*d + (1:d);
            H(idx, idx) = H_i;
        end
    end
end
