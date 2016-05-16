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
blocklength = 300000;
for i=1:length(SNRs)    
    [~, pe_le(i), ~] = estimate_pbit(@simulate_le, SNRs(i), err_needed, blocklength);
        
    [~, pe_dfe(i), ~] = estimate_pbit(@simulate_dfe, SNRs(i), err_needed, blocklength);
    
    [~, pe_aa(i), ~] = estimate_pbit(@simulate_aa, SNRs(i), err_needed, blocklength);
        
    [~, pe_aa2(i), ~] = estimate_pbit(@simulate_aa2, SNRs(i), err_needed, blocklength);
    
    [~, pe_sim(i), ~] = estimate_pbit(@simulation_bound, SNRs(i), err_needed, blocklength);
    
%    [~, pe_vit(i), ~] = estimate_pbit(@simulate_viterbi, SNRs(i), err_needed, blocklength);
    
%     fprintf('.');
%     if mod(i, 5) == 0
%         fprintf('%d', i);
%     end
end
pes(:, 1) = pe_le;
pes(:, 2) = pe_dfe;
pes(:, 3) = pe_aa;
pes(:, 4) = pe_aa2;
pes(:, 5) = pe_vit;
pes(:, 7) = pe_sim;

% fprintf('\n');

% Plots
figure;
semilogy(SNRs, pes(:,1), 'b:');
hold on;
semilogy(SNRs, pes(:,2), 'b-');
semilogy(SNRs, pes(:,3), 'r:');
semilogy(SNRs, pes(:,4), 'r-');
semilogy(SNRs, pes(:,5), 'c:');
semilogy(SNRs, pes(:,7), 'c');

semilogy(SNRs, pe_bounds(:,1), 'g-');

xlabel('SNR [dB]');
ylabel('Pe');
grid on;
xlim([8, 14]);
ylim([1e-4, 1e-1]);
