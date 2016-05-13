function [pbit, pe, tot_bits] = estimate_pbit(simulate_txrx, SNR, err_needed, blocklength)
tot_errs = 0;
tot_sym_errs = 0;
tot_bits = 0;
while tot_errs < err_needed
    [~, ~, err_count, sym_err_count] = simulate_txrx(blocklength, SNR);
    tot_bits = tot_bits + blocklength;
    tot_errs = tot_errs + err_count;
    tot_sym_errs = tot_sym_errs + sym_err_count;
end
pbit = tot_errs / tot_bits;
pe = tot_sym_errs / (tot_bits/2);
end
