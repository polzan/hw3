close all; clear all;

alphabet = [1+1j; 1-1j; -1-1j; -1+1j];
%alphabet = [];
M = length(alphabet);
K = 20;
L1 = 1;
L2 = 2;

states = all_states(alphabet, L1+L2);

path_metrics = sparse(M^(L1+L2), K+1);
survivors = sparse(M^(L1+L2), K+1);

% Initial state set to -1,-1
%path_metrics(:,1) = [0; Inf * ones(M^(L1+L2)-1,1)];
% First state is the initial state k=-1
rng(4);
rho = alphabet(round(rand(K, 1) .* 3 + 1));
f = @(sigma_j, sigma_i) sigma_j(1+L1);
for k=0:K-1
    next_path_metrics = zeros(length(states), 1);
    next_survivors = zeros(length(states), 1);
   for j=1:length(states)
       best_i = 0;
       best_path_metric = Inf;
       for i=1:length(states)
           previous_path_metric = full(path_metrics(i, (k-1)+2));
           if previous_path_metric == Inf; continue; end
           branch_metr_i_j = branch_metric(states(j,:), states(i,:), rho(k+1), f);
           if branch_metr_i_j == Inf; continue; end
           next_path_metric = previous_path_metric + branch_metr_i_j;
           if next_path_metric < best_path_metric
               best_i = i;
               best_path_metric = next_path_metric;
           end
       end
       next_path_metrics(j) = best_path_metric;
       next_survivors(j) = best_i;
   end
   path_metrics(:, k+2) = next_path_metrics;
   survivors(:, k+2) = next_survivors;
end

full(path_metrics)
full(survivors)

% Detected state sequence
[min_metric, min_j] = min(full(path_metrics(:, K+1)));
state_seq = zeros(K+1, L1+L2);
state_seq(K-1+2,:) = states(min_j,:);
k = K-1;
j = min_j;
while k > -1
    previous_state_i = full(survivors(j, k+2));
    previous_state = states(previous_state_i,:);
    state_seq((k-1)+2,:) = previous_state;
    k = k-1;
    j = previous_state_i;
end

% Detected symbols
a = zeros(K, 1);
k = K-1;
while k >= 0
    a(k+1) = state_seq(k+2, L1+1);
    k = k-1;
end

sym_err_count = sum(a ~= rho);
Pe = sym_err_count / K;

fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);

