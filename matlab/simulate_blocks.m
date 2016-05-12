close all; clear all; clc;

SNRs = linspace(8, 14, 20);

pbits = zeros(length(SNRs), 2);
pbit_bounds = zeros(length(SNRs), 1);
for i=1:length(SNRs)    
    [pbits(i, 1), tot_errs, tot_bits] = estimate_pbit(@simulate_le, SNRs(i));
    fprintf('LE Nbit = %d, errs = %d\n', tot_bits, tot_errs);
    
    [pbits(i, 2), tot_errs, tot_bits] = estimate_pbit(@simulate_dfe, SNRs(i));
    fprintf('DFE Nbit = %d, errs = %d\n', tot_bits, tot_errs);
    
    SNR_lin = 10^(SNRs(i)/10);
    Pe = 4*(1-1/sqrt(4))*(1 - normcdf(sqrt(3/(4-1)*SNR_lin), 0, 1));
    pbit_bounds(i,1) = 1/log2(4) * Pe;
end

figure;

semilogy(SNRs, pbits(:,1));
hold on;
semilogy(SNRs, pbits(:,2));
semilogy(SNRs, pbit_bounds(:,1));
