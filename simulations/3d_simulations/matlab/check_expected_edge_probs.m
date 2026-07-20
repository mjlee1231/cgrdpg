% Check what edge probabilities we SHOULD expect with NEW design
clear; clc;

fprintf('Expected Edge Probabilities with NEW Design\n');
fprintf('============================================\n\n');

n = 1000;
tau = 0.001;

% NEW latent design
t = (1:n)' / n;
X0 = [0.3*t + 0.5, 0.15 * sin(2*pi*t) + 0.6, 0.1 * cos(4*pi*t)];
S = diag([1, 1, -1]);

fprintf('Latent position ranges:\n');
fprintf('  Dim 1: [%.4f, %.4f]\n', min(X0(:,1)), max(X0(:,1)));
fprintf('  Dim 2: [%.4f, %.4f]\n', min(X0(:,2)), max(X0(:,2)));
fprintf('  Dim 3: [%.4f, %.4f]\n\n', min(X0(:,3)), max(X0(:,3)));

% Compute all pairwise inner products
Y0 = X0 * S;  % Y = X * S
LinearPred = X0 * Y0';  % (n × d) * (d × n) = (n × n)

% Extract upper triangle (actual edges)
triu_idx = triu(true(n), 1);
s_values = LinearPred(triu_idx);

fprintf('Linear predictor s_ij = x_i^T * S * x_j:\n');
fprintf('  Min:    %.4f\n', min(s_values));
fprintf('  Median: %.4f\n', median(s_values));
fprintf('  Max:    %.4f\n\n', max(s_values));

% Edge probabilities
eps_clip = 1e-10;
p_values = max(min(1 ./ (1 + exp(-s_values / tau)), 1 - eps_clip), eps_clip);

fprintf('Edge probabilities with tau = %.4f:\n', tau);
fprintf('  Min:    %.10f\n', min(p_values));
fprintf('  Median: %.10f\n', median(p_values));
fprintf('  Max:    %.10f\n', max(p_values));
fprintf('  # exactly 1.0 (after clip): %d / %d\n', sum(p_values >= 1 - eps_clip), length(p_values));

% Check a few specific values
fprintf('\nSample edge probabilities:\n');
for k = 1:5
    i = k;
    j = k + 1;
    s_ij = X0(i,:) * S * X0(j,:)';
    p_ij = 1 / (1 + exp(-s_ij / tau));
    fprintf('  p(%d,%d): s=%.4f, p=%.10f\n', i, j, s_ij, p_ij);
end

fprintf('\nConclusion:\n');
if min(p_values) > 0.99
    fprintf('  ⚠️  ALL edge probs > 0.99 with tau=%.4f\n', tau);
    fprintf('  This might be TOO saturated for stable optimization\n');
    fprintf('  Consider using larger tau (e.g., 0.01 or 0.1)\n');
else
    fprintf('  Edge probabilities have reasonable range\n');
end
