close all;

clc;

alphabet = [-1,1];

K = 30;
L1 = 1;
L2 = 1;

psi = [0; -0.2; 1; -0.2; 0];
D = 2;
rng(5);
a_tx = round(rand(1e4, 1)) .* 2 - 1;
w = sqrt(2/10^(11/10)) .* (randn(length(a_tx), 1) + 1j .* randn(length(a_tx), 1));

rho = filter(psi, 1, a_tx)+w;

detector = ViterbiDetector(alphabet, K, L1, L2, psi, D);

a = detector.detect(rho);

output_translation = K-1;

a_cut_first = a((output_translation+D+1:length(a)));
sym_err_count = sum(a_cut_first ~= a_tx((1:length(a_cut_first))));
Pe = sym_err_count / length(a_cut_first);

fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);

% figure;
% stem(0:length(rho)-1, real(rho));
% hold on;
% stem(0:length(a)-1, real(a));
% ylim([-1.2 1.2]);

figure;
% subplot(1,2,1);
stem((0:length(a_tx)-1) + D +output_translation, 1.1.*real(a_tx));
hold on;
stem(0:length(a)-1, real(a));
ylim([-1.2 1.2]);

% subplot(1,2,2);
% stem(0:length(a_tx)-1, 1.1.*imag(a_tx));
% hold on;
% stem(0:length(a)-1, imag(a));
% ylim([-1.2 1.2]);

% figure;
% scatter(real(rho), imag(rho));
% hold on;
% scatter(real(a), imag(a));
