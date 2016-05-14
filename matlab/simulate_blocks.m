close all; clear all; clc;

SNRs = linspace(8, 14, 20);

% Symbol error bounds
pe_bounds = zeros(length(SNRs), 2);
for i=1:length(SNRs)
    Q = @(x) 1 - normcdf(x, 0, 1);
    
    % QAM lower bound
    SNR_lin = 10^(SNRs(i)/10);    
    pe_bounds(i,1) = 2 * Q(sqrt(SNR_lin));
end

% Pe vs SNR
pes = zeros(length(SNRs), 2);
err_needed = 30;
blocklength = 1e4;
for i=1:length(SNRs)    
    [~, pes(i, 1), tot_bits] = estimate_pbit(@simulate_le, SNRs(i), err_needed, blocklength);
    %fprintf('LE Nbit = %d\n', tot_bits);
    
    [~, pes(i, 2), tot_bits] = estimate_pbit(@simulate_dfe, SNRs(i), err_needed, blocklength);
    %fprintf('DFE Nbit = %d, errs = %d\n', tot_bits, tot_errs);
    
    [~, pes(i, 3), tot_bits] = estimate_pbit(@simulate_viterbi, SNRs(i), err_needed, blocklength);
    
    fprintf('.');
    if mod(i, 5) == 0
        fprintf('%d', i);
    end
end
fprintf('\n');

% Plots
figure;
semilogy(SNRs, pes(:,1), 'b:');
hold on;
semilogy(SNRs, pes(:,2), 'b-');
semilogy(SNRs, pes(:,3), 'r:');

semilogy(SNRs, pe_bounds(:,1), 'g-');

xlabel('SNR [dB]');
ylabel('Pe');
grid on;
xlim([8, 14]);
ylim([1e-4, 1e-1]);
