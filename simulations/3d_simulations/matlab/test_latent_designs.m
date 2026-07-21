% Quick test of different latent position designs
clear; clc;

fprintf('Testing Latent Position Designs\n');
fprintf('================================\n\n');

n = 1000;
S = diag([1, 1, -1]);

designs = {'current', 'scaled', 'orthogonal'};
design_names = {'Current (original)', 'Scaled (2x coefficients)', 'Orthogonalized'};

fprintf('%-25s %8s %8s %8s %10s %10s %8s\n', ...
    'Design', 'λ1', 'λ2', 'λ3', 'Gap(1-2)', 'Gap(2-3)', 'Cond#');
fprintf('%s\n', repmat('-', 1, 85));

for i = 1:length(designs)
    X0 = generate_latent_positions(n, designs{i});

    % Eigenvalues
    G = X0' * X0;
    eigs = sort(eig(G), 'descend');
    gap1 = eigs(1) - eigs(2);
    gap2 = eigs(2) - eigs(3);
    cond_num = cond(G);

    % Edge probabilities
    P = X0 * S * X0';
    p_min = min(P(:));
    p_max = max(P(:));

    fprintf('%-25s %8.2f %8.2f %8.2f %10.2f %10.2f %8.2f\n', ...
        design_names{i}, eigs(1), eigs(2), eigs(3), gap1, gap2, cond_num);
    fprintf('  Edge prob range: [%.4f, %.4f]\n', p_min, p_max);
end

fprintf('\n');
fprintf('Criteria check:\n');
fprintf('  ✓ Good eigenvalue gap: > 50\n');
fprintf('  ✓ Large smallest eigenvalue: > 20\n');
fprintf('  ✓ Good condition number: < 20\n');
fprintf('  ✓ Valid edge probabilities: [0, 1]\n');
