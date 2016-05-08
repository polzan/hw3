close all; clear all; clc;

T = 1;
T_Q = T/4;
qc_length = 32;
Qf_samples = 1000;

[b,a] = transmitter_tf();
q_c = impz(b, a, qc_length);
[Q_f_half, f_half] = freqz(b, a, Qf_samples, 1/T_Q);

figure;
stem((0:length(q_c)-1), q_c);
xlabel('k');
ylabel('q_c(k)');

figure;
plot(f_half, 20*log10(abs(Q_f_half)));
xlabel('f (Fs = 4/T)');
ylabel('|Q_c(f)| [dB]');
