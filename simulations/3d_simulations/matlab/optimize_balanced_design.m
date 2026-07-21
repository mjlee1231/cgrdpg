% Optimize balanced design to keep edge probabilities in [0,1]
% while maintaining good eigenvalue properties
clear; clc;

fprintf('Optimizing Balanced Design for Valid Edge Probabilities\n');
fprintf('========================================================\n\n');

n = 1000;
t = (1:n)' / n;
S = diag([1, 1, -1]);

% Starting point: Balanced design that exceeded 1.0
% [0.4*t + 0.55, 0.25*sin(2πt) + 0.55, 0.18*cos(4πt)]

% Try different scaling and shift combinations
% Format: [scale_t, shift_t, scale_sin, shift_sin, scale_cos]
candidates = [
    % Reduce shifts
    0.40, 0.50, 0.25, 0.50, 0.18;
    0.40, 0.48, 0.25, 0.48, 0.18;
    0.40, 0.45, 0.25, 0.45, 0.18;

    % Reduce oscillation scales slightly
    0.38, 0.50, 0.23, 0.50, 0.16;
    0.36, 0.52, 0.22, 0.52, 0.15;

    % More conservative
    0.35, 0.52, 0.20, 0.52, 0.14;
    0.38, 0.48, 0.24, 0.48, 0.17;

    % Try to maximize λ3 while keeping P valid
    0.40, 0.47, 0.26, 0.47, 0.19;
    0.42, 0.46, 0.27, 0.46, 0.20;
];

fprintf('Testing %d candidate designs...\n\n', size(candidates, 1));

results = cell(size(candidates, 1), 1);

for i = 1:size(candidates, 1)
    scale_t = candidates(i, 1);
    shift_t = candidates(i, 2);
    scale_sin = candidates(i, 3);
    shift_sin = candidates(i, 4);
    scale_cos = candidates(i, 5);

    X0 = [scale_t * t + shift_t, ...
          scale_sin * sin(2*pi*t) + shift_sin, ...
          scale_cos * cos(4*pi*t)];

    % Compute eigenvalues
    G = X0' * X0;
    eigs = sort(eig(G), 'descend');
    gap1 = eigs(1) - eigs(2);
    gap2 = eigs(2) - eigs(3);
    cond_num = cond(G);

    % Compute edge probabilities
    P = X0 * S * X0';
    p_min = min(P(:));
    p_max = max(P(:));

    % Check criteria
    valid_p = (p_min >= 0) && (p_max <= 1.0);
    good_gap2 = gap2 >= 20;
    good_lambda3 = eigs(3) >= 10;
    good_cond = cond_num <= 100;

    all_criteria = valid_p && good_gap2 && good_lambda3 && good_cond;

    results{i} = struct('idx', i, 'params', candidates(i, :), ...
        'eigs', eigs, 'gap1', gap1, 'gap2', gap2, 'cond', cond_num, ...
        'p_min', p_min, 'p_max', p_max, ...
        'valid_p', valid_p, 'good_gap2', good_gap2, ...
        'good_lambda3', good_lambda3, 'good_cond', good_cond, ...
        'all_criteria', all_criteria);
end

% Display results
fprintf('%-4s %-6s %-6s %-6s %-6s %-6s | %-8s %-8s %-8s | %-8s %-8s | %-8s %-8s | %s\n', ...
    'Idx', 'sc_t', 'sh_t', 'sc_si', 'sh_si', 'sc_co', ...
    'λ1', 'λ2', 'λ3', 'gap2', 'cond#', 'P_min', 'P_max', 'Valid?');
fprintf('%s\n', repmat('-', 1, 130));

for i = 1:length(results)
    r = results{i};
    valid_str = 'X';
    if r.all_criteria
        valid_str = 'OK';
    end
    fprintf('%-4d %.2f   %.2f   %.2f   %.2f   %.2f  | %8.2f %8.2f %8.2f | %8.2f %8.2f | %8.4f %8.4f | %s\n', ...
        r.idx, r.params(1), r.params(2), r.params(3), r.params(4), r.params(5), ...
        r.eigs(1), r.eigs(2), r.eigs(3), r.gap2, r.cond, ...
        r.p_min, r.p_max, valid_str);
end

% Find best valid design
valid_designs = cellfun(@(x) x.all_criteria, results);
if any(valid_designs)
    fprintf('\n========================================\n');
    fprintf('VALID DESIGNS (meet all criteria):\n');
    fprintf('========================================\n\n');

    valid_idx = find(valid_designs);

    % Sort by λ3 (largest smallest eigenvalue)
    lambda3_vals = cellfun(@(x) x.eigs(3), results(valid_idx));
    [~, sort_idx] = sort(lambda3_vals, 'descend');

    for k = 1:length(valid_idx)
        i = valid_idx(sort_idx(k));
        r = results{i};

        fprintf('Design %d (Rank %d by λ3):\n', r.idx, k);
        fprintf('  X0 = [%.2f*t + %.2f, %.2f*sin(2πt) + %.2f, %.2f*cos(4πt)]\n', ...
            r.params(1), r.params(2), r.params(3), r.params(4), r.params(5));
        fprintf('  Eigenvalues: λ1=%.2f, λ2=%.2f, λ3=%.2f\n', r.eigs(1), r.eigs(2), r.eigs(3));
        fprintf('  Gaps: λ1-λ2=%.2f, λ2-λ3=%.2f\n', r.gap1, r.gap2);
        fprintf('  Condition number: %.2f\n', r.cond);
        fprintf('  Edge prob range: [%.4f, %.4f]\n\n', r.p_min, r.p_max);
    end

    % Recommend the best
    best_idx = valid_idx(sort_idx(1));
    best = results{best_idx};

    fprintf('========================================\n');
    fprintf('RECOMMENDED DESIGN (Design %d):\n', best.idx);
    fprintf('========================================\n');
    fprintf('X0 = [%.2f*t + %.2f, %.2f*sin(2πt) + %.2f, %.2f*cos(4πt)]\n\n', ...
        best.params(1), best.params(2), best.params(3), best.params(4), best.params(5));
    fprintf('Summary:\n');
    fprintf('  * λ3 = %.2f (largest among valid designs)\n', best.eigs(3));
    fprintf('  * λ2 - λ3 = %.2f > 20\n', best.gap2);
    fprintf('  * Condition# = %.2f < 100\n', best.cond);
    fprintf('  * Edge prob range: [%.4f, %.4f] in [0, 1]\n', best.p_min, best.p_max);

else
    fprintf('\n! No designs meet all criteria. Showing best partial matches:\n\n');

    % Find designs with valid edge probabilities
    valid_p_designs = cellfun(@(x) x.valid_p, results);
    if any(valid_p_designs)
        fprintf('Designs with valid edge probabilities:\n');
        for i = find(valid_p_designs)'
            r = results{i};
            fprintf('  Design %d: λ3=%.2f, gap2=%.2f, cond=%.2f, P=[%.4f, %.4f]\n', ...
                r.idx, r.eigs(3), r.gap2, r.cond, r.p_min, r.p_max);
        end
    end
end
