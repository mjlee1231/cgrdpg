% Design better latent positions with well-separated eigenvalues
% Keep edge probability range similar to current [0.487942, 0.997058]
% Test different parameter combinations

clear; clc;

fprintf('Testing Different Latent Position Designs\n');
fprintf('==========================================\n\n');

n = 1000;
t = (1:n)' / n;
S = diag([1, 1, -1]);

% Target edge probability range
target_min = 0.48;
target_max = 0.99;

fprintf('Target: Edge prob range [%.3f, %.3f]\n', target_min, target_max);
fprintf('Goal: Better eigenvalue separation\n\n');

%% Define candidate designs
designs = {};

% Current design
designs{1} = struct(...
    'name', 'Current (3D helix)', ...
    'formula', '0.15*sin(2πt)+0.6, 0.15*cos(2πt)+0.6, 0.15*cos(4πt)', ...
    'X', [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t)]);

% Design 2: Larger amplitude in dim 1
designs{2} = struct(...
    'name', 'Larger amp dim1', ...
    'formula', '0.25*sin(2πt)+0.6, 0.1*cos(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.25 * sin(2*pi*t) + 0.6, ...
          0.1 * cos(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

% Design 3: Different frequencies
designs{3} = struct(...
    'name', 'Freq 1,2,3', ...
    'formula', '0.2*sin(2πt)+0.6, 0.15*cos(4πt)+0.6, 0.1*cos(6πt)', ...
    'X', [0.2 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(4*pi*t) + 0.6, ...
          0.1 * cos(6*pi*t)]);

% Design 4: Unequal amplitudes with offset differences
designs{4} = struct(...
    'name', 'Unequal amps+offsets', ...
    'formula', '0.2*sin(2πt)+0.7, 0.12*cos(2πt)+0.55, 0.08*cos(4πt)', ...
    'X', [0.2 * sin(2*pi*t) + 0.7, ...
          0.12 * cos(2*pi*t) + 0.55, ...
          0.08 * cos(4*pi*t)]);

% Design 5: Linear + periodic
designs{5} = struct(...
    'name', 'Linear+periodic', ...
    'formula', '0.3*t+0.5, 0.15*sin(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.3*t + 0.5, ...
          0.15 * sin(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

% Design 6: Scaled dim 3 only
designs{6} = struct(...
    'name', 'Scaled dim3', ...
    'formula', '0.15*sin(2πt)+0.6, 0.15*cos(2πt)+0.6, 0.25*cos(4πt)', ...
    'X', [0.15 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.25 * cos(4*pi*t)]);

% Design 7: All different amplitudes
designs{7} = struct(...
    'name', 'Amps 0.2,0.15,0.1', ...
    'formula', '0.2*sin(2πt)+0.6, 0.15*cos(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.2 * sin(2*pi*t) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

% Design 8: Step-like in dim 1
designs{8} = struct(...
    'name', 'Step+periodic', ...
    'formula', '0.2*sign(sin(2πt))+0.6, 0.15*cos(2πt)+0.6, 0.1*cos(4πt)', ...
    'X', [0.2 * sign(sin(2*pi*t)) + 0.6, ...
          0.15 * cos(2*pi*t) + 0.6, ...
          0.1 * cos(4*pi*t)]);

%% Evaluate each design
fprintf('%-3s | %-25s | %-15s | %-25s | %-10s\n', ...
    '#', 'Design', 'Edge Range', 'Top 3 Eigenvalues', 'Ratios');
fprintf('%s\n', repmat('-', 1, 100));

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

    % Ratios
    ratio_21 = abs(top3(2)) / abs(top3(1));
    ratio_32 = abs(top3(3)) / abs(top3(2));

    % Check if in target range
    in_range = (p_min >= target_min - 0.05) && (p_max <= target_max + 0.05);

    % Print
    fprintf('%-3d | %-25s | [%.3f,%.3f]', ...
        i, designs{i}.name, p_min, p_max);

    if in_range
        fprintf(' ✓ ');
    else
        fprintf(' ✗ ');
    end

    fprintf(' | %+7.2f,%+7.2f,%+7.2f | %.3f,%.3f', ...
        top3(1), top3(2), top3(3), ratio_21, ratio_32);

    fprintf('\n');

    % Store results
    results(i).design_num = i;
    results(i).name = designs{i}.name;
    results(i).p_min = p_min;
    results(i).p_max = p_max;
    results(i).in_range = in_range;
    results(i).eigvals = top3;
    results(i).ratio_21 = ratio_21;
    results(i).ratio_32 = ratio_32;
    results(i).X = X;
end

fprintf('%s\n', repmat('-', 1, 100));

%% Find best designs
fprintf('\n========================================\n');
fprintf('Best Designs (in target range)\n');
fprintf('========================================\n\n');

% Filter designs in range
valid_designs = find([results.in_range]);

if isempty(valid_designs)
    fprintf('⚠ No designs in target edge prob range!\n');
    fprintf('Consider adjusting parameters.\n');
else
    % Sort by eigenvalue separation (want smaller ratios)
    separation_scores = [results(valid_designs).ratio_21] + [results(valid_designs).ratio_32];
    [~, sort_idx] = sort(separation_scores);

    fprintf('Ranked by eigenvalue separation (lower ratio = better):\n\n');

    for rank = 1:min(3, length(valid_designs))
        idx = valid_designs(sort_idx(rank));
        r = results(idx);

        fprintf('%d. %s\n', rank, r.name);
        fprintf('   Edge range: [%.4f, %.4f]\n', r.p_min, r.p_max);
        fprintf('   Eigenvalues: %+.3f, %+.3f, %+.3f\n', r.eigvals);
        fprintf('   Ratios: λ₂/λ₁=%.3f, λ₃/λ₂=%.3f (sum=%.3f)\n', ...
            r.ratio_21, r.ratio_32, r.ratio_21 + r.ratio_32);
        fprintf('   Formula: %s\n\n', designs{idx}.formula);
    end

    % Recommend best
    best_idx = valid_designs(sort_idx(1));
    fprintf('========================================\n');
    fprintf('RECOMMENDED DESIGN: #%d - %s\n', best_idx, results(best_idx).name);
    fprintf('========================================\n\n');

    fprintf('MATLAB code:\n');
    fprintf('X0 = %s;\n\n', designs{best_idx}.formula);

    % Show improvement
    current_idx = 1;
    if best_idx ~= current_idx
        fprintf('Improvement over current:\n');
        fprintf('  Current ratios: λ₂/λ₁=%.3f, λ₃/λ₂=%.3f (sum=%.3f)\n', ...
            results(current_idx).ratio_21, results(current_idx).ratio_32, ...
            results(current_idx).ratio_21 + results(current_idx).ratio_32);
        fprintf('  New ratios:     λ₂/λ₁=%.3f, λ₃/λ₂=%.3f (sum=%.3f)\n', ...
            results(best_idx).ratio_21, results(best_idx).ratio_32, ...
            results(best_idx).ratio_21 + results(best_idx).ratio_32);

        improvement = (results(current_idx).ratio_21 + results(current_idx).ratio_32) - ...
                     (results(best_idx).ratio_21 + results(best_idx).ratio_32);
        fprintf('  Separation improvement: %.3f (%.1f%%)\n', improvement, ...
            100 * improvement / (results(current_idx).ratio_21 + results(current_idx).ratio_32));
    end
end

%% Detailed analysis of top 3
fprintf('\n========================================\n');
fprintf('Detailed Analysis of Top 3 Candidates\n');
fprintf('========================================\n\n');

for rank = 1:min(3, length(valid_designs))
    idx = valid_designs(sort_idx(rank));
    X = designs{idx}.X;
    Y = X * S;
    P = X * Y';

    fprintf('Design #%d: %s\n', idx, designs{idx}.name);
    fprintf('%s\n', repmat('-', 1, 50));

    % Eigenvalue analysis
    eigvals_all = eig(P);
    [~, sort_idx_all] = sort(abs(eigvals_all), 'descend');
    top_10 = eigvals_all(sort_idx_all(1:min(10, length(eigvals_all))));

    fprintf('Top 10 eigenvalues:\n');
    for j = 1:length(top_10)
        fprintf('  λ_%d = %+10.4f', j, top_10(j));
        if j <= 3
            fprintf(' ← embedding dimension');
        end
        fprintf('\n');
    end

    % Spectral gaps
    fprintf('\nSpectral gaps:\n');
    fprintf('  Gap 1-2: %.4f\n', abs(top_10(1)) - abs(top_10(2)));
    fprintf('  Gap 2-3: %.4f\n', abs(top_10(2)) - abs(top_10(3)));
    fprintf('  Gap 3-4: %.4f\n', abs(top_10(3)) - abs(top_10(4)));

    % Condition number
    fprintf('\nNumerical properties:\n');
    fprintf('  Condition number of P:     %.2e\n', cond(P));
    fprintf('  Condition number of X^T X: %.2e\n', cond(X'*X));

    fprintf('\n');
end
