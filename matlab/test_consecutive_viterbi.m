close all;

alphabet = [1+1j; 1-1j; -1-1j; -1+1j];
%alphabet = [-1,1];

K = 20;
L1 = 0;
L2 = 2;

psi = [0; 0; 1; 0.1];
t0 = 1;
D = 1;

rng(4);
a_tx = alphabet(round(rand(3*K, 1) .* 3 + 1));

rho = conv(psi, a_tx);

initial_path_metrics = zeros(4^(L1+L2), 1);
a_total = [];
for i=1:3
    rho_block = rho(K*(i-1)+(1:K));
    [a, final_path_metrics] = viterbi_decoder(L1, L2, K, psi, t0, D, initial_path_metrics, rho_block);
    
    % Check error on this block
    a_cut_first = a(t0+D+1:length(a));
    sym_err_count = sum(a_cut_first ~= a_tx(K*(i-1)+(1:length(a_cut_first))));
    Pe = sym_err_count / K;
    fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);
    
    a_total = [a_total; a];
    intial_path_metrics = final_path_metrics;
end

% Check total errors
a_tot_cut_first = a_total(t0+D+1:length(a_total));
sym_err_count = sum(a_tot_cut_first ~= a_tx((1:length(a_tot_cut_first))));
Pe = sym_err_count / K;
fprintf('Total Pe = %f\n', Pe);

figure;
hold on;
stem(0:length(rho)-1, real(rho));
stem(0:length(a_total)-1, real(a_total));


