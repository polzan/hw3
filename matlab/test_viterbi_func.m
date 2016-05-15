close all; 

alphabet = [1+1j; 1-1j; -1-1j; -1+1j];
%alphabet = [-1,1];

K = 300;
L1 = 0;
L2 = 2;

psi = [0; 0; 1.4; 0; 0];
t0 = 0;
D = 2;

initial_path_metrics = zeros(4^(L1+L2), 1);

rng(4);
a_tx = alphabet(round(rand(K, 1) .* 3 + 1));
w = sqrt(0.1) .* (randn(length(a_tx), 1) + 1j .* randn(length(a_tx), 1));

rho = conv(psi, a_tx + w);

detector = ViterbiDetector(alphabet, K, L1, L2, psi, t0 + D);


[a, final_path_metrics] = detector.detect_symbols(rho, initial_path_metrics);

output_translation = L2+L1;

a_cut_first = a(max(output_translation+1,1):min(length(a),length(a)+output_translation));

sym_err_count = sum(a_cut_first ~= a_tx(1:length(a_cut_first)));
Pe = sym_err_count / K;

fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);

figure;
stem(0:length(rho)-1, real(rho));
hold on;
stem(0:length(a)-1, real(a));
ylim([-1.2 1.2]);

figure;
subplot(1,2,1);
stem(0:length(a_tx)-1, 1.1.*real(a_tx));
hold on;
stem(0:length(a)-1, real(a));
ylim([-1.2 1.2]);

subplot(1,2,2);
stem(0:length(a_tx)-1, 1.1.*imag(a_tx));
hold on;
stem(0:length(a)-1, imag(a));
ylim([-1.2 1.2]);

figure;
scatter(real(rho), imag(rho));
hold on;
scatter(real(a), imag(a));


