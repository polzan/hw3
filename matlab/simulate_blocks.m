close all; clear all; clc;

SNRs = linspace(8, 14, 5); % odd number so it uses SNR=11

% Symbol error bounds
pe_bounds = zeros(length(SNRs), 2);
for i=1:length(SNRs)
    Q = @(x) 1 - normcdf(x, 0, 1);
    
    % QAM lower bound
    SNR_lin = 10^(SNRs(i)/10);
    pe_bounds(i,1) = 2 * Q(sqrt(SNR_lin));
end

% Pe vs SNR
pes = zeros(length(SNRs), 8);
err_needed = 30;
blocklength = 1e5;
for i=1:length(SNRs)
    [~, pe_le(i), ~] = estimate_pbit(@simulate_le, SNRs(i), err_needed, blocklength);
    
    [~, pe_dfe(i), ~] = estimate_pbit(@simulate_dfe, SNRs(i), err_needed, blocklength);
    
    [~, pe_aa(i), ~] = estimate_pbit(@simulate_aa, SNRs(i), err_needed, blocklength);
    
    [~, pe_aa2(i), ~] = estimate_pbit(@simulate_aa2, SNRs(i), err_needed, blocklength);
    
    [~, pe_sim(i), ~] = estimate_pbit(@simulation_bound, SNRs(i), err_needed, blocklength);
    
    [~, pe_vit(i), ~] = estimate_pbit(@simulate_viterbi, SNRs(i), err_needed, blocklength);
    
    [~, pe_mlm(i), ~] = estimate_pbit(@simulate_maxlogmap, SNRs(i), err_needed, blocklength);
    
    fprintf('.');
    if mod(i, 5) == 0
        fprintf('%d', i);
    end
end

fprintf('\n');

% Plots
figure;
semilogy(SNRs, pe_le, 'b:');
hold on;
semilogy(SNRs, pe_dfe, 'b-');
semilogy(SNRs, pe_aa2, 'k:');
semilogy(SNRs, pe_aa, 'k-');
semilogy(SNRs, pe_vit, 'r:');
semilogy(SNRs, pe_mlm, 'r');

semilogy(SNRs, pe_bounds(:,1), 'g-');
semilogy(SNRs, pe_sim, 'g:');

xlabel('SNR [dB]');
ylabel('Pe');
grid on;
xlim([8, 14]);
ylim([1e-4, 1e-1]);


save('pe_plots_workspace.mat');
