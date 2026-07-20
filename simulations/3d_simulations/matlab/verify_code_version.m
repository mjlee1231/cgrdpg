% Verify the code has the bug fix applied
clear; clc;

fprintf('Code Version Check\n');
fprintf('==================\n\n');

% Check if the fix is present by examining the function
type core/compute_cgrdpg_one_step_batch.m

fprintf('\n\n');
fprintf('Look for these lines around line 33-40:\n');
fprintf('  Y_i = S_estimated * x_i;  (CORRECT)\n');
fprintf('  s_i = X_init(idx_j, :) * Y_i;  (CORRECT)\n\n');
fprintf('If you see instead:\n');
fprintf('  Y_init = X_init * S_estimated;  (WRONG)\n');
fprintf('  s_i = Y_init(idx_j, :) * x_i;  (WRONG)\n\n');
fprintf('Then the code was not updated!\n');
