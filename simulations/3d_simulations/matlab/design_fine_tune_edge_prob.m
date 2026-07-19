% Fine-tune latent positions to achieve max edge prob ~0.995-0.997
% while maintaining well-separated eigenvalues
%
% Strategy: Start from current design and carefully increase offsets/amplitudes
%           to push max edge prob from 0.990 → 0.997 without exceeding 1.0

clear; clc;

fprintf('Fine-tuning Latent Positions for Max Edge Prob ~0.997\n');
fprintf('======================================================\n\n');

n = 1000;
t = (1:n)' / n;
S = diag([1, 1, -1]);

% Target edge probability range
target_max = 0.997;
tolerance = 0.005;  % Accept 0.992-0.998

fprintf('Target: Max edge prob = %.3f ± %.3f\n', target_max, tolerance);
fprintf('Constraint: Edge prob must be ≤ 1.0\n');
fprintf('Goal: Eigenvalue gap > 5\n\n');

%% Define candidate designs
% Start from current and make SMALL adjustments
designs = {};

% Reference: Current design
designs{1} = struct(...
    'name', 'CURRENT', ...
    'formula', '0.3*t+0.5, 0.15*sin(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.5, ...
          0.15 * sin(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

% Design 2: Increase linear offset slightly
designs{2} = struct(...
    'name', 'Linear +0.02 offset', ...
    'formula', '0.3*t+0.52, 0.15*sin(2πt)+0.62, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.52, ...
          0.15 * sin(2*pi*t) + 0.62, ...
          0.1 * cos(4*pi*t)]);

% Design 3: Increase linear offset more
designs{3} = struct(...
    'name', 'Linear +0.03 offset', ...
    'formula', '0.3*t+0.53, 0.15*sin(2πt)+0.63, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.53, ...
          0.15 * sin(2*pi*t) + 0.63, ...
          0.1 * cos(4*pi*t)]);

% Design 4: Increase periodic amplitude
designs{4} = struct(...
    'name', 'Periodic +0.02 amp', ...
    'formula', '0.3*t+0.52, 0.17*sin(2πt)+0.62, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.52, ...
          0.17 * sin(2*pi*t) + 0.62, ...
          0.1 * cos(4*pi*t)]);

% Design 5: Increase linear slope
designs{5} = struct(...
    'name', 'Linear slope 0.32', ...
    'formula', '0.32*t+0.52, 0.15*sin(2πt)+0.62, 0.1*cos(4πt)', ...
    'X', [0.32*t + 0.52, ...
          0.15 * sin(2*pi*t) + 0.62, ...
          0.1 * cos(4*pi*t)]);

% Design 6: Increase linear slope more
designs{6} = struct(...
    'name', 'Linear slope 0.33', ...
    'formula', '0.33*t+0.52, 0.15*sin(2πt)+0.62, 0.1*cos(4πt)', ...
    'X', [0.33*t + 0.52, ...
          0.15 * sin(2*pi*t) + 0.62, ...
          0.1 * cos(4*pi*t)]);

% Design 7: Both slope and offset increase
designs{7} = struct(...
    'name', 'Slope 0.32 + offset 0.53', ...
    'formula', '0.32*t+0.53, 0.16*sin(2πt)+0.62, 0.1*cos(4πt)', ...
    'X', [0.32*t + 0.53, ...
          0.16 * sin(2*pi*t) + 0.62, ...
          0.1 * cos(4*pi*t)]);

% Design 8: Larger increase
designs{8} = struct(...
    'name', 'Slope 0.33 + offset 0.54', ...
    'formula', '0.33*t+0.54, 0.16*sin(2πt)+0.63, 0.1*cos(4πt)', ...
    'X', [0.33*t + 0.54, ...
          0.16 * sin(2*pi*t) + 0.63, ...
          0.1 * cos(4*pi*t)]);

% Design 9: Aggressive increase (testing limit)
designs{9} = struct(...
    'name', 'Slope 0.35 + offset 0.54', ...
    'formula', '0.35*t+0.54, 0.17*sin(2πt)+0.63, 0.09*cos(4πt)', ...
    'X', [0.35*t + 0.54, ...
          0.17 * sin(2*pi*t) + 0.63, ...
          0.09 * cos(4*pi*t)]);

% Design 10: OLD design for comparison
designs{10} = struct(...
    'name', 'OLD (3D helix)', ...
    'formula', '0.15*sin(2πt)+0.6, 0.15*cos(2πt)+0.6, 0.15*cos(4πt)', ...
    'X', [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)]);

% Design 11: Modified helix with better separation
designs{11} = struct(...
    'name', 'Modified helix (unequal amps)', ...
    'formula', '0.2*sin(2πt)+0.6, 0.13*cos(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.2 * sin(2*pi*t) + 0.6, ...
          0.13 * cos(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

% Design 12: Even more unequal for separation
designs{12} = struct(...
    'name', 'Helix very unequal amps', ...
    'formula', '0.22*sin(2πt)+0.6, 0.12*cos(2πt)+0.6, 0.08*cos(4πt)', ...
    'X', [0.22 * sin(2*pi*t) + 0.6, ...
          0.12 * cos(2*pi*t) + 0.6, ...
          0.08 * cos(4*pi*t)]);

%% Evaluate each design
fprintf('%-3s | %-30s | %-18s | %-30s | %-10s | %-8s\n', ...
    '#', 'Design', 'Edge Range', 'Top 3 Eigenvalues', 'Gap 2-3', 'Valid?');
fprintf('%s\n', repmat('-', 1, 125));

results = [];

for i = 1:length(designs)
    X = designs{i}.X;
    Y = X * S;
    P = X * Y';

    % Edge prob range
    p_min = min(P(:));
    p_max = max(P(:));

    % Eigenvalues
    eigvals = eig(P);
    [~, idx] = sort(abs(eigvals), 'descend');
    top3 = eigvals(idx(1:3));

    % Gap
    gap_23 = abs(abs(top3(2)) - abs(top3(3)));

    % Check validity
    valid_range = (p_max >= target_max - tolerance) && (p_max <= target_max + tolerance);
    valid_prob = (p_min >= 0) && (p_max <= 1.0);
    well_separated = gap_23 > 5;

    % Print
    fprintf('%-3d | %-30s | [%.4f,%.4f]', i, designs{i}.name, p_min, p_max);

    if valid_range
        fprintf(' ✓ ');
    else
        fprintf('   ');
    end

    fprintf(' | %+7.2f,%+7.2f,%+7.2f | %.2f', top3(1), top3(2), top3(3), gap_23);

    if well_separated
        fprintf(' ✓   ');
    else
        fprintf('     ');
    end

    fprintf('| ');
    if valid_prob && well_separated
        fprintf('✓✓');
    elseif valid_prob
        fprintf('✓ ');
    else
        fprintf('✗ (>1)');
    end

    fprintf('\n');

    % Store results
    results(i).design_num = i;
    results(i).name = designs{i}.name;
    results(i).formula = designs{i}.formula;
    results(i).p_min = p_min;
    results(i).p_max = p_max;
    results(i).valid_range = valid_range;
    results(i).valid_prob = valid_prob;
    results(i).eigvals = top3;
    results(i).gap_23 = gap_23;
    results(i).well_separated = well_separated;
    results(i).X = X;
    results(i).score = 0;

    % Scoring: how close to target max edge prob
    if valid_prob && well_separated
        results(i).score = -abs(p_max - target_max);  % Negative distance (higher = better)
    end
end

fprintf('%s\n', repmat('-', 1, 125));

%% Find best designs
fprintf('\n========================================\n');
fprintf('Best Designs\n');
fprintf('========================================\n\n');

% Filter: valid probability range AND well separated
good_designs = find([results.valid_prob] & [results.well_separated]);

if isempty(good_designs)
    fprintf('⚠ No valid designs found!\n');
else
    % Sort by score (closest to target max edge prob)
    scores = [results(good_designs).score];
    [~, sort_idx] = sort(scores, 'descend');

    fprintf('Ranked by proximity to target max edge prob (%.3f):\n\n', target_max);

    for rank = 1:min(5, length(good_designs))
        idx = good_designs(sort_idx(rank));
        r = results(idx);

        fprintf('%d. %s\n', rank, r.name);
        fprintf('   Edge range: [%.4f, %.4f]\n', r.p_min, r.p_max);
        fprintf('   Distance to target: %.4f\n', abs(r.p_max - target_max));
        fprintf('   Eigenvalues: %+.3f, %+.3f, %+.3f\n', r.eigvals);
        fprintf('   Gap 2-3: %.2f\n', r.gap_23);
        fprintf('   Formula: %s\n\n', r.formula);
    end

    % Recommend best
    best_idx = good_designs(sort_idx(1));
    fprintf('========================================\n');
    fprintf('RECOMMENDED: Design #%d - %s\n', best_idx, results(best_idx).name);
    fprintf('========================================\n\n');

    fprintf('Edge prob range: [%.4f, %.4f]\n', ...
        results(best_idx).p_min, results(best_idx).p_max);
    fprintf('Eigenvalues: %+.3f, %+.3f, %+.3f\n', results(best_idx).eigvals);
    fprintf('Gap 2-3: %.2f\n\n', results(best_idx).gap_23);

    fprintf('MATLAB code:\n');
    fprintf('X0 = [%s];\n\n', results(best_idx).formula);

    % Compare with current (#1) and old (#10)
    current_idx = 1;
    old_idx = 10;

    fprintf('Comparison:\n');
    fprintf('  Old max edge prob:     %.4f (gap=%.2f)\n', ...
        results(old_idx).p_max, results(old_idx).gap_23);
    fprintf('  Current max edge prob: %.4f (gap=%.2f)\n', ...
        results(current_idx).p_max, results(current_idx).gap_23);
    fprintf('  NEW max edge prob:     %.4f (gap=%.2f)\n', ...
        results(best_idx).p_max, results(best_idx).gap_23);

    fprintf('\nImprovement vs current:\n');
    fprintf('  Max edge prob: +%.4f (+%.1f%%)\n', ...
        results(best_idx).p_max - results(current_idx).p_max, ...
        100 * (results(best_idx).p_max - results(current_idx).p_max) / results(current_idx).p_max);

    if results(best_idx).gap_23 > results(old_idx).gap_23
        fprintf('  Eigenvalue separation: MUCH BETTER than old (gap %.2f vs %.2f)\n', ...
            results(best_idx).gap_23, results(old_idx).gap_23);
    end
end
