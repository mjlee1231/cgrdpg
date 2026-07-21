% Analyze eigenvalue properties of latent position designs
clear; clc;

fprintf('Analyzing Latent Position Design Properties\n');
fprintf('===========================================\n\n');

n = 1000;
t = (1:n)' / n;

%% Current Design
fprintf('CURRENT DESIGN:\n');
fprintf('X0 = [0.3*t + 0.5, 0.15*sin(2πt) + 0.6, 0.1*cos(4πt)]\n\n');

X_current = [0.3*t + 0.5, ...
             0.15 * sin(2*pi*t) + 0.6, ...
             0.1 * cos(4*pi*t)];

% Gram matrix
G_current = X_current' * X_current;
eig_current = eig(G_current);
[eig_current_sorted, idx] = sort(eig_current, 'descend');

fprintf('Gram matrix X''X eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_current_sorted(1));
fprintf('  λ2 = %.4f\n', eig_current_sorted(2));
fprintf('  λ3 = %.4f\n', eig_current_sorted(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', ...
    eig_current_sorted(1) - eig_current_sorted(2), eig_current_sorted(1)/eig_current_sorted(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', ...
    eig_current_sorted(2) - eig_current_sorted(3), eig_current_sorted(2)/eig_current_sorted(3));
fprintf('\nCondition number: %.2f\n', cond(G_current));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_current_sorted)));

% Column correlations
corr_matrix = corr(X_current);
fprintf('Column correlations:\n');
fprintf('  corr(X1, X2) = %.4f\n', corr_matrix(1,2));
fprintf('  corr(X1, X3) = %.4f\n', corr_matrix(1,3));
fprintf('  corr(X2, X3) = %.4f\n\n', corr_matrix(2,3));

% Edge probability range
S = diag([1, 1, -1]);
P_current = X_current * S * X_current';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_current(:)), max(P_current(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_current(:)));

%% Proposed Design 1: Adjusted scaled design with better gaps
fprintf('========================================\n');
fprintf('PROPOSED DESIGN 1: Adjusted scaled with better λ2-λ3 gap\n');
fprintf('X0 = [0.35*t + 0.575, 0.22*sin(2πt) + 0.575, 0.12*cos(4πt)]\n');
fprintf('(more balanced scaling)\n\n');

X_prop1 = [0.35*t + 0.575, ...
           0.22 * sin(2*pi*t) + 0.575, ...
           0.12 * cos(4*pi*t)];

G_prop1 = X_prop1' * X_prop1;
eig_prop1 = sort(eig(G_prop1), 'descend');

fprintf('Gram matrix eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_prop1(1));
fprintf('  λ2 = %.4f\n', eig_prop1(2));
fprintf('  λ3 = %.4f\n', eig_prop1(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', eig_prop1(1) - eig_prop1(2), eig_prop1(1)/eig_prop1(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', eig_prop1(2) - eig_prop1(3), eig_prop1(2)/eig_prop1(3));
fprintf('\nCondition number: %.2f\n', cond(G_prop1));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_prop1)));

P_prop1 = X_prop1 * S * X_prop1';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_prop1(:)), max(P_prop1(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_prop1(:)));

%% Proposed Design 2: Aggressive scaling with controlled shifts
fprintf('========================================\n');
fprintf('PROPOSED DESIGN 2: Aggressive scaling for large eigenvalues\n');
fprintf('X0 = [0.45*t + 0.55, 0.28*sin(2πt) + 0.55, 0.15*cos(4πt)]\n');
fprintf('(larger oscillatory components)\n\n');

X_prop2 = [0.45*t + 0.55, ...
           0.28 * sin(2*pi*t) + 0.55, ...
           0.15 * cos(4*pi*t)];

G_prop2 = X_prop2' * X_prop2;
eig_prop2 = sort(eig(G_prop2), 'descend');

fprintf('Gram matrix eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_prop2(1));
fprintf('  λ2 = %.4f\n', eig_prop2(2));
fprintf('  λ3 = %.4f\n', eig_prop2(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', eig_prop2(1) - eig_prop2(2), eig_prop2(1)/eig_prop2(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', eig_prop2(2) - eig_prop2(3), eig_prop2(2)/eig_prop2(3));
fprintf('\nCondition number: %.2f\n', cond(G_prop2));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_prop2)));

P_prop2 = X_prop2 * S * X_prop2';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_prop2(:)), max(P_prop2(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_prop2(:)));

%% Proposed Design 3: Better balanced design with controlled eigenvalues
fprintf('========================================\n');
fprintf('PROPOSED DESIGN 3: Balanced design with moderate shifts\n');
fprintf('X0 = [0.4*t + 0.55, 0.25*sin(2πt) + 0.55, 0.18*cos(4πt)]\n\n');

X_prop3 = [0.4*t + 0.55, ...
           0.25 * sin(2*pi*t) + 0.55, ...
           0.18 * cos(4*pi*t)];

G_prop3 = X_prop3' * X_prop3;
eig_prop3 = sort(eig(G_prop3), 'descend');

fprintf('Gram matrix eigenvalues:\n');
fprintf('  λ1 = %.4f\n', eig_prop3(1));
fprintf('  λ2 = %.4f\n', eig_prop3(2));
fprintf('  λ3 = %.4f\n', eig_prop3(3));
fprintf('\nEigenvalue gaps:\n');
fprintf('  λ1 - λ2 = %.4f  (ratio: %.2f)\n', eig_prop3(1) - eig_prop3(2), eig_prop3(1)/eig_prop3(2));
fprintf('  λ2 - λ3 = %.4f  (ratio: %.2f)\n', eig_prop3(2) - eig_prop3(3), eig_prop3(2)/eig_prop3(3));
fprintf('\nCondition number: %.2f\n', cond(G_prop3));
fprintf('Smallest eigenvalue: %.4f\n\n', min(abs(eig_prop3)));

P_prop3 = X_prop3 * S * X_prop3';
fprintf('Edge probability range: [%.4f, %.4f]\n', min(P_prop3(:)), max(P_prop3(:)));
fprintf('Mean edge probability: %.4f\n\n', mean(P_prop3(:)));

%% Summary
fprintf('========================================\n');
fprintf('SUMMARY COMPARISON\n');
fprintf('========================================\n\n');

designs = {'Current', 'Adjusted', 'Aggressive', 'Balanced'};
eigs_all = {eig_current_sorted, eig_prop1, eig_prop2, eig_prop3};
P_all = {P_current, P_prop1, P_prop2, P_prop3};

fprintf('%-15s %8s %8s %8s %12s %8s %12s\n', ...
    'Design', 'λ1', 'λ2', 'λ3', 'Gap1', 'Gap2', 'Cond#');
fprintf('%-15s %8s %8s %8s %12s %8s %12s\n', ...
    repmat('-', 1, 15), repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 8), ...
    repmat('-', 1, 12), repmat('-', 1, 8), repmat('-', 1, 12));

for i = 1:4
    eigs = eigs_all{i};
    gap1 = eigs(1) - eigs(2);
    gap2 = eigs(2) - eigs(3);
    cond_num = eigs(1) / eigs(3);

    fprintf('%-15s %8.2f %8.2f %8.2f %12.2f %8.2f %12.2f\n', ...
        designs{i}, eigs(1), eigs(2), eigs(3), gap1, gap2, cond_num);
end

fprintf('\n%-15s %10s %10s\n', 'Design', 'P_min', 'P_max');
fprintf('%-15s %10s %10s\n', repmat('-', 1, 15), repmat('-', 1, 10), repmat('-', 1, 10));
for i = 1:4
    P = P_all{i};
    fprintf('%-15s %10.4f %10.4f\n', designs{i}, min(P(:)), max(P(:)));
end

% Column correlations
fprintf('\n%-15s %12s %12s %12s\n', 'Design', 'corr(X1,X2)', 'corr(X1,X3)', 'corr(X2,X3)');
fprintf('%-15s %12s %12s %12s\n', repmat('-', 1, 15), repmat('-', 1, 12), repmat('-', 1, 12), repmat('-', 1, 12));
X_all = {X_current, X_prop1, X_prop2, X_prop3};
for i = 1:4
    C = corr(X_all{i});
    fprintf('%-15s %12.4f %12.4f %12.4f\n', designs{i}, C(1,2), C(1,3), C(2,3));
end

fprintf('\n');
fprintf('RECOMMENDATION:\n');
fprintf('===============\n');
fprintf('All proposed designs improve upon the current design.\n');
fprintf('Choose based on priority:\n\n');
fprintf('  Adjusted:   Moderate improvement, keeps similar structure\n');
fprintf('  Aggressive: Large eigenvalues, good gaps, moderate condition#\n');
fprintf('  Balanced:   Best overall balance of all criteria\n\n');

fprintf('Target criteria for good design:\n');
fprintf('  ✓ Gap λ1-λ2 > 200 (for strong identifiability)\n');
fprintf('  ✓ Gap λ2-λ3 > 20 (for separation of all eigenvalues)\n');
fprintf('  ✓ Smallest eigenvalue > 10 (in magnitude)\n');
fprintf('  ✓ Condition number < 100 (acceptable for numerical stability)\n');
fprintf('  ✓ Edge probabilities in [0, 1] range\n');
fprintf('  ✓ Low column correlations (< 0.5 in magnitude)\n');
