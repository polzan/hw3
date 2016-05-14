close all; 

alphabet = [1+1j; 1-1j; -1-1j; -1+1j];
%alphabet = [-1,1];

K = 100;
L1 = 0;
L2 = 2;

psi = [0; 0; 1; 0.1];
t0 = 1;
D = 1;

initial_path_metrics = zeros(4^(L1+L2), 1);

rng(4);
a_tx = alphabet(round(rand(K, 1) .* 3 + 1));

rho = conv(psi, a_tx);
%rho = rho(t0+D+1:length(rho)-1);

[a, final_path_metrics] = viterbi_decoder(L1, L2, K, psi, t0, D, initial_path_metrics, rho);

a_cut_first = a(t0+D+1:length(a));

sym_err_count = sum(a_cut_first ~= a_tx(1:length(a_cut_first)));
Pe = sym_err_count / K;

fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);

figure;
stem(0:length(rho)-1, real(rho));
hold on;
stem(0:length(a)-1, real(a));
ylim([-1.2 1.2]);

figure;
scatter(real(rho), imag(rho));
hold on;
scatter(real(a), imag(a));
