function [pbit, tot_errs, tot_bits] = estimate_pbit(simulate_txrx, SNR)
err_needed = 100;
blocklength = 1e4;
tot_errs = 0;
tot_bits = 0;
while tot_errs < err_needed
    [~, err_count] = simulate_txrx(blocklength, SNR);
    tot_bits = tot_bits + blocklength;
    tot_errs = tot_errs + err_count;
end
pbit = tot_errs / tot_bits;
end
