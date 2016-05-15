close all;

clc;

alphabet = [-1,1];

K = 50;
L1 = 1;
L2 = 3;

psi = [0.75; 0.01];
D = 0;
rng(5);
a_tx = round(rand(1e3, 1)) .* 2 - 1;

rho = filter(psi, 1, a_tx);

detector = ViterbiDetector(alphabet, K, L1, L2, psi, D);

a = detector.detect(rho);

output_translation = K-1;

a_cut_first = a(output_translation+1:length(a));

sym_err_count = sum(a_cut_first ~= a_tx(1:length(a_cut_first)));
Pe = sym_err_count / K;

fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);

% figure;
% stem(0:length(rho)-1, real(rho));
% hold on;
% stem(0:length(a)-1, real(a));
% ylim([-1.2 1.2]);

figure;
% subplot(1,2,1);
stem((0:length(a_tx)-1) + output_translation, 1.1.*real(a_tx));
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
