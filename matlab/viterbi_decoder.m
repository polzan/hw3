function [detected_syms, final_path_metrics] = viterbi_decoder(L1, L2, trellis_depth, psi, t0, D, initial_path_metrics, rho)
alphabet = [1+1j; 1-1j; -1-1j; -1+1j]; % QPSK

M = length(alphabet);
K = trellis_depth;
total_psi_delay = t0+D;

states = all_states(alphabet, L1+L2);

path_metrics = sparse(M^(L1+L2), K+1);
survivors = sparse(M^(L1+L2), K+1);

path_metrics(:,1) = initial_path_metrics; % First state is the initial state k=-1

f = @(sigma_j, sigma_i) received_sample(sigma_j, sigma_i, psi, total_psi_delay, L1, L2);

for k=0:K-1
    next_path_metrics = zeros(length(states), 1);
    next_survivors = zeros(length(states), 1);
    for j=1:length(states)
        best_i = 0;
        best_path_metric = Inf;
        for i=1:length(states)
            previous_path_metric = full(path_metrics(i, (k-1)+2));
            if previous_path_metric == Inf; continue; end
            branch_metr_i_j = branch_metric(j, i, states(j,:), states(i,:), rho(k+1), f);
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

% Detected state sequence
final_path_metrics = full(path_metrics(:, K+1));
[~, min_j] = min(final_path_metrics);
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
detected_syms = zeros(K, 1);
k = K-1;
while k >= 0
    detected_syms(k+1) = state_seq(k+2, L1+1);
    k = k-1;
end
end

function s = all_states(alphabet, L)
if L > 1
    s_ = all_states(alphabet, L-1);
    s = [];
    for i=1:length(alphabet)
        for j=1:size(s_,1)
            r = [s_(j,:), alphabet(i)];
            s = [s; r];
        end
    end
elseif L == 1
    for i=1:length(alphabet)
        s(i,1) = alphabet(i);
    end
else
    error('L<1');
end
end

function bm = branch_metric(j, i, sigma_j, sigma_i, rho_k, f)
if ~is_transition_possible(j,i,sigma_j,sigma_i)
    bm = Inf;
else
    bm = abs(rho_k - f(sigma_j, sigma_i))^2;
end
end

function y = is_transition_possible(j, i, sigma_j, sigma_i)
persistent result_table;
if isempty(result_table)
    result_table = sparse(1000, 1000);
end
y_ = result_table(j, i);
if y_ == 0
    y_ = 1 + any(sigma_j(2:length(sigma_j)) ~= sigma_i(1:length(sigma_i)-1));
    result_table(j, i) = y_;
end
y = y_ - 1;
end

function u_k = received_sample(sigma_j, sigma_i, psi, psi_delay, L1, L2)
a = [sigma_j(1,:), sigma_i(1,length(sigma_i))];
assert(length(a) == L1 + L2 + 1);
u = conv(a, psi);
u_k = u(psi_delay+L1+1);
end
