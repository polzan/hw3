close all; clear all; clc;

T = 1;
T_Q = T/4;
Qf_samples = 1000;

[b,a,qc_length,q_c] = transmitter_tf();
[Q_f_half, f_half] = freqz(b, a, Qf_samples, 1/T_Q);

figure;
stem((0:length(q_c)-1), q_c);
xlabel('n [@ T/4]');
ylabel('q_c(nT/4)');
grid on;
print('plot_qc', '-depsc');

figure;
plot(f_half, 20*log10(abs(Q_f_half)));
xlabel('f (@ 4/T)');
ylabel('|Q_c(f)| [dB]');
grid on;
print('plot_Qf', '-depsc');