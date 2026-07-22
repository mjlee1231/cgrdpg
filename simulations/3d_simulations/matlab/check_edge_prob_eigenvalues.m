% Check eigenvalues of edge probability matrix for optimized design
clear; clc;

fprintf('Checking eigenvalues of P = X0 * S * X0'' for optimized design\n');
fprintf('===============================================================\n\n');

n = 1000;
t = (1:n)' / n;
S = diag([1, 1, -1]);

% Optimized design
X0 = [0.42 * t + 0.46, ...
      0.27 * sin(2*pi*t) + 0.46, ...
      0.20 * cos(4*pi*t)];

fprintf('Design: X0 = [0.42*t + 0.46, 0.27*sin(2πt) + 0.46, 0.20*cos(4πt)]\n');
fprintf('Signature: S = diag(1, 1, -1)\n\n');

% Gram matrix (what we optimized)
G = X0' * X0;
eigs_G = sort(eig(G), 'descend');
fprintf('Eigenvalues of Gram matrix G = X0''*X0:\n');
fprintf('  λ1(G) = %.2f\n', eigs_G(1));
fprintf('  λ2(G) = %.2f\n', eigs_G(2));
fprintf('  λ3(G) = %.2f\n\n', eigs_G(3));

% Edge probability matrix
P = X0 * S * X0';
eigs_P = sort(eig(P), 'descend');

% Find non-zero eigenvalues (rank is at most 3)
nonzero_eigs = eigs_P(abs(eigs_P) > 1e-10);

fprintf('Eigenvalues of edge probability matrix P = X0*S*X0'':\n');
fprintf('  Number of non-zero eigenvalues: %d\n', length(nonzero_eigs));
fprintf('  Non-zero eigenvalues:\n');
for i = 1:length(nonzero_eigs)
    fprintf('    λ%d(P) = %+.2f\n', i, nonzero_eigs(i));
end

fprintf('\n');
fprintf('Expected signature: (+, +, -) from S = diag(1, 1, -1)\n');
fprintf('Actual signs: (%s, %s, %s)\n', ...
    sign_str(nonzero_eigs(1)), ...
    sign_str(nonzero_eigs(2)), ...
    sign_str(nonzero_eigs(3)));

% Edge probability range
fprintf('\nEdge probability range:\n');
fprintf('  min(P) = %.4f\n', min(P(:)));
fprintf('  max(P) = %.4f\n', max(P(:)));

% Check if P has correct indefinite structure
fprintf('\nVerification:\n');
if length(nonzero_eigs) == 3 && nonzero_eigs(1) > 0 && nonzero_eigs(2) > 0 && nonzero_eigs(3) < 0
    fprintf('  ✓ P has correct indefinite structure: 2 positive, 1 negative eigenvalue\n');
else
    fprintf('  ✗ WARNING: P does not have expected indefinite structure!\n');
end

% Alternative: Check eigenvalues of the 3x3 matrix X0'*P*X0 = X0'*X0*S*X0'*X0
fprintf('\nAlternative check - eigenvalues of G*S*G (should match signature):\n');
GSG = G * S * G;
eigs_GSG = sort(eig(GSG), 'descend');
for i = 1:3
    fprintf('  λ%d(GSG) = %+.2f\n', i, eigs_GSG(i));
end

function s = sign_str(x)
    if x > 0
        s = '+';
    elseif x < 0
        s = '-';
    else
        s = '0';
    end
end
