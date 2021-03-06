close all; clear all; clc;

load('pe_plots_data.mat');

% Plots
figure;
semilogy(SNRs, pe_le, 'b--');
hold on;
semilogy(SNRs, pe_dfe, 'b-');
semilogy(SNRs, pe_aa2, 'k--');
semilogy(SNRs, pe_aa, 'k-');
semilogy(SNRs, pe_vit, 'r--');
semilogy(SNRs, pe_mlm, 'r');

semilogy(SNRs, pe_noisi, 'g--');
semilogy(SNRs, pe_theor_bound, 'g-');

xlabel('SNR [dB]');
ylabel('Pe');
grid on;
xlim([8, 14]);
ylim([1e-4, 1e-1]);

legend('MF+LE@T', 'MF+DFE@T', 'AAF+DFE@T/2', 'AAF+DFE@T', 'VA', 'FBA', 'MFb-s', 'MFb-t');

print('pe_plot', '-depsc');
