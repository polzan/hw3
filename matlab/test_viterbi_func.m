close all; 

alphabet = [1+1j; 1-1j; -1-1j; -1+1j];
%alphabet = [-1,1];

K = 100;
L1 = 0;
L2 = 3;

psi = [0; 0; 1; 0.1];
t0 = 1;
D = 1;

initial_path_metrics = zeros(4^(L1+L2), 1);

rng(4);
a_tx = alphabet(round(rand(K, 1) .* 3 + 1));

rho = conv(psi, a_tx);

detector = ViterbiDetector(alphabet, K, L1, L2, psi, t0 + D);


[a, final_path_metrics] = detector.detect_symbols(rho, initial_path_metrics);

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
