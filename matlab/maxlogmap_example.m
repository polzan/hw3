close all;

clc;

alphabet = [-1,1];

K = 20;
Kin = 600;
L1 = 0;
L2 = 2;

psi = 1;
D = 0;
rng(5);
a_tx = round(rand(1e4, 1)) .* 2 - 1;
w = sqrt(2/10^(11/10)) .* (randn(length(a_tx), 1) + 1j .* randn(length(a_tx), 1));

rho = filter(psi, 1, a_tx)+w;


detector = MaxLogMapDetector(alphabet, K, Kin, L1, L2, psi, D);

a = detector.detect(rho);

drop_output = K-1;

a_cut = a(drop_output+1:length(a)-drop_output);
a_tx_cut = a_tx((drop_output+1:length(a_tx)-drop_output)-D);
sym_err_count = sum(a_cut ~= a_tx_cut);
Pe = sym_err_count / length(a_cut);

fprintf('Pe = %f on one trellis depth of %d\n', Pe, K);

% figure;
% stem(0:length(rho)-1, real(rho));
% hold on;
% stem(0:length(a)-1, real(a));
% ylim([-1.2 1.2]);

figure;
% subplot(1,2,1);
stem(0:length(a_tx)-1, 0.9.*real(a_tx));
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
