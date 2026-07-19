% Design latent positions with HIGH edge probabilities (max ~0.997)
% while maintaining WELL-SEPARATED eigenvalues
%
% Current issue:
%   Old design: max edge prob = 0.9971, but degenerate eigenvalues (gap=0)
%   New design: eigenvalues separated (gap=11.66), but max edge prob = 0.99 only
%
% Goal: Combine best of both - high edge prob + separated eigenvalues

clear; clc;

fprintf('Designing Latent Positions: High Edge Prob + Separated Eigenvalues\n');
fprintf('===================================================================\n\n');

n = 1000;
t = (1:n)' / n;
S = diag([1, 1, -1]);

% Target edge probability range
target_min = 0.48;
target_max_low = 0.990;   % Minimum acceptable max
target_max_high = 0.998;  % Maximum acceptable max

fprintf('Target: Max edge prob in [%.3f, %.3f]\n', target_max_low, target_max_high);
fprintf('Goal: Eigenvalue gap > 5 (well-separated)\n\n');

%% Define candidate designs
% Strategy: Increase offsets to push edge probs higher
designs = {};

% Reference: Old design (degenerate but high edge prob)
designs{1} = struct(...
    'name', 'OLD: 3D helix (degenerate)', ...
    'formula', '0.15*sin(2πt)+0.6, 0.15*cos(2πt)+0.6, 0.15*cos(4πt)', ...
    'X', [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)]);

% Reference: Current design (separated but low edge prob)
designs{2} = struct(...
    'name', 'CURRENT: Linear+periodic', ...
    'formula', '0.3*t+0.5, 0.15*sin(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.5, ...
          0.15 * sin(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

% Design 3: Linear with higher offset
designs{3} = struct(...
    'name', 'Linear high offset', ...
    'formula', '0.3*t+0.55, 0.15*sin(2πt)+0.65, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.55, ...
          0.15 * sin(2*pi*t) + 0.65, ...
          0.1 * cos(4*pi*t)]);

% Design 4: Linear with even higher offset
designs{4} = struct(...
    'name', 'Linear higher offset', ...
    'formula', '0.3*t+0.6, 0.15*sin(2πt)+0.7, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.6, ...
          0.15 * sin(2*pi*t) + 0.7, ...
          0.1 * cos(4*pi*t)]);

% Design 5: Larger linear slope + high offset
designs{5} = struct(...
    'name', 'Large slope + high offset', ...
    'formula', '0.4*t+0.5, 0.15*sin(2πt)+0.65, 0.08*cos(4πt)', ...
    'X', [0.4*t + 0.5, ...
          0.15 * sin(2*pi*t) + 0.65, ...
          0.08 * cos(4*pi*t)]);

% Design 6: Periodic with higher offsets
designs{6} = struct(...
    'name', 'Periodic high offsets', ...
    'formula', '0.2*sin(2πt)+0.65, 0.15*cos(2πt)+0.7, 0.1*cos(4πt)+0.05', ...
    'X', [0.2 * sin(2*pi*t) + 0.65, ...
          0.15 * cos(2*pi*t) + 0.7, ...
          0.1 * cos(4*pi*t) + 0.05]);

% Design 7: Unequal amplitudes + very high offsets
designs{7} = struct(...
    'name', 'Unequal amps + very high offsets', ...
    'formula', '0.25*sin(2πt)+0.7, 0.12*cos(2πt)+0.65, 0.08*cos(4πt)', ...
    'X', [0.25 * sin(2*pi*t) + 0.7, ...
          0.12 * cos(2*pi*t) + 0.65, ...
          0.08 * cos(4*pi*t)]);

% Design 8: Linear quadratic blend
designs{8} = struct(...
    'name', 'Linear-quadratic blend', ...
    'formula', '0.3*t+0.2*t.^2+0.4, 0.15*sin(2πt)+0.65, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.2*t.^2 + 0.4, ...
          0.15 * sin(2*pi*t) + 0.65, ...
          0.1 * cos(4*pi*t)]);

% Design 9: Different frequency + high offsets
designs{9} = struct(...
    'name', 'Freq variation + high offsets', ...
    'formula', '0.25*sin(2πt)+0.65, 0.15*cos(4πt)+0.7, 0.1*cos(6πt)', ...
    'X', [0.25 * sin(2*pi*t) + 0.65, ...
          0.15 * cos(4*pi*t) + 0.7, ...
          0.1 * cos(6*pi*t)]);

% Design 10: Push to max ~0.997 with linear
designs{10} = struct(...
    'name', 'Push max to 0.997 (linear)', ...
    'formula', '0.35*t+0.6, 0.15*sin(2πt)+0.7, 0.08*cos(4πt)', ...
    'X', [0.35*t + 0.6, ...
          0.15 * sin(2*pi*t) + 0.7, ...
          0.08 * cos(4*pi*t)]);

%% Evaluate each design
fprintf('%-3s | %-30s | %-18s | %-30s | %-15s\n', ...
    '#', 'Design', 'Edge Range', 'Top 3 Eigenvalues', 'Gap 2-3');
fprintf('%s\n', repmat('-', 1, 120));

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

    % Ratios and gap
    ratio_21 = abs(top3(2)) / abs(top3(1));
    ratio_32 = abs(top3(3)) / abs(top3(2));
    gap_23 = abs(abs(top3(2)) - abs(top3(3)));

    % Check if in target range
    in_max_range = (p_max >= target_max_low) && (p_max <= target_max_high);
    well_separated = gap_23 > 5;

    % Print
    fprintf('%-3d | %-30s | [%.4f,%.4f]', i, designs{i}.name, p_min, p_max);

    if in_max_range
        fprintf(' ✓ ');
    else
        fprintf(' ✗ ');
    end

    fprintf(' | %+7.2f,%+7.2f,%+7.2f | %.2f', top3(1), top3(2), top3(3), gap_23);

    if well_separated
        fprintf(' ✓');
    else
        fprintf(' ✗');
    end

    fprintf('\n');

    % Store results
    results(i).design_num = i;
    results(i).name = designs{i}.name;
    results(i).formula = designs{i}.formula;
    results(i).p_min = p_min;
    results(i).p_max = p_max;
    results(i).in_max_range = in_max_range;
    results(i).eigvals = top3;
    results(i).ratio_21 = ratio_21;
    results(i).ratio_32 = ratio_32;
    results(i).gap_23 = gap_23;
    results(i).well_separated = well_separated;
    results(i).X = X;
end

fprintf('%s\n', repmat('-', 1, 120));

%% Find best designs
fprintf('\n========================================\n');
fprintf('Best Designs (high edge prob + separated eigenvalues)\n');
fprintf('========================================\n\n');

% Filter: in max range AND well separated
good_designs = find([results.in_max_range] & [results.well_separated]);

if isempty(good_designs)
    fprintf('⚠ No designs satisfy both criteria!\n');
    fprintf('Relaxing constraints...\n\n');

    % Show designs with good separation (even if max edge prob not perfect)
    separated_designs = find([results.well_separated]);
    if ~isempty(separated_designs)
        fprintf('Designs with good eigenvalue separation (gap > 5):\n\n');
        for idx = separated_designs
            r = results(idx);
            fprintf('%d. %s\n', idx, r.name);
            fprintf('   Edge range: [%.4f, %.4f]\n', r.p_min, r.p_max);
            fprintf('   Eigenvalues: %+.3f, %+.3f, %+.3f\n', r.eigvals);
            fprintf('   Gap 2-3: %.2f\n', r.gap_23);
            fprintf('   Formula: %s\n\n', r.formula);
        end
    end

    % Show designs with high max edge prob (even if not well separated)
    high_edge_designs = find([results.in_max_range]);
    if ~isempty(high_edge_designs)
        fprintf('Designs with high max edge prob (%.3f-%.3f):\n\n', ...
            target_max_low, target_max_high);
        for idx = high_edge_designs
            r = results(idx);
            fprintf('%d. %s\n', idx, r.name);
            fprintf('   Edge range: [%.4f, %.4f]\n', r.p_min, r.p_max);
            fprintf('   Eigenvalues: %+.3f, %+.3f, %+.3f\n', r.eigvals);
            fprintf('   Gap 2-3: %.2f\n', r.gap_23);
            fprintf('   Formula: %s\n\n', r.formula);
        end
    end
else
    % Sort by gap (larger gap = better separation)
    gaps = [results(good_designs).gap_23];
    [~, sort_idx] = sort(gaps, 'descend');

    fprintf('Ranked by eigenvalue gap (higher = better separation):\n\n');

    for rank = 1:length(good_designs)
        idx = good_designs(sort_idx(rank));
        r = results(idx);

        fprintf('%d. %s\n', rank, r.name);
        fprintf('   Edge range: [%.4f, %.4f]\n', r.p_min, r.p_max);
        fprintf('   Eigenvalues: %+.3f, %+.3f, %+.3f\n', r.eigvals);
        fprintf('   Gap 2-3: %.2f  (ratio λ₃/λ₂ = %.3f)\n', r.gap_23, r.ratio_32);
        fprintf('   Formula: %s\n\n', r.formula);
    end

    % Recommend best
    best_idx = good_designs(sort_idx(1));
    fprintf('========================================\n');
    fprintf('RECOMMENDED: Design #%d - %s\n', best_idx, results(best_idx).name);
    fprintf('========================================\n\n');

    fprintf('MATLAB code:\n');
    fprintf('X0 = %s;\n\n', results(best_idx).formula);

    % Compare with current
    current_idx = 2;  % Design 2 is current
    fprintf('Improvement over current design:\n');
    fprintf('  Current max edge prob: %.4f\n', results(current_idx).p_max);
    fprintf('  New max edge prob:     %.4f (%.2f%% higher)\n', ...
        results(best_idx).p_max, ...
        100 * (results(best_idx).p_max - results(current_idx).p_max) / results(current_idx).p_max);
    fprintf('  Current gap 2-3:       %.2f\n', results(current_idx).gap_23);
    fprintf('  New gap 2-3:           %.2f\n', results(best_idx).gap_23);
end
