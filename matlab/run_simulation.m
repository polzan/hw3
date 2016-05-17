close all; clear all; clc;

Nsyms = 3e5;
Nbits = 2*Nsyms;

% Use more points where they are useful
SNRs = linspace(8, 12, 20);
SNRs = [SNRs, 13, 14];
% Always compute for SNR=11
if isempty(find(SNRs == 11, 1))
    warning('Add SNR = 11 dB. The plots will have one more point');
    SNRs = [SNRs(SNRs < 11), 11, SNRs(SNRs > 11)];
end

% QAM lower bound
Q = @(x) 1 - normcdf(x, 0, 1);
SNRs_lin = power(10, SNRs ./ 10);
pe_theor_bound = 2 .* Q(sqrt(SNRs_lin));

% Pe vs SNR
for i=1:length(SNRs)
    [pbit_le(i), pe_le(i)] = simulate_le(Nbits, SNRs(i));
    [pbit_dfe(i), pe_dfe(i)] = simulate_dfe(Nbits, SNRs(i));
    [pbit_aa(i), pe_aa(i)] = simulate_aa(Nbits, SNRs(i));
    [pbit_aa2(i), pe_aa2(i)] = simulate_aa2(Nbits, SNRs(i));
    [pbit_vit(i), pe_vit(i)] = simulate_viterbi(Nbits, SNRs(i));
    [pbit_mlm(i), pe_mlm(i)] = simulate_maxlogmap(Nbits, SNRs(i));
    
    [pbit_noisi(i), pe_noisi(i)] = simulation_bound(Nbits, SNRs(i));
    
    if mod(i, 5) == 0
        fprintf('%d', i);
    else
        fprintf('.');
    end
end
fprintf('\n');
save('pe_plots_data.mat');
