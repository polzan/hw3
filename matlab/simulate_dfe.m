function [Pbit, Pe, err_count, sym_err_count] = simulate_dfe(Nbits, SNR)
% Parameters
t0 = 33;
t0_sampled = floor(t0/4);
M1 = 3;
M2 = 3;
D = 1;

Nbits_w_transient = Nbits + 2*t0_sampled + 2*D;

[rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(Nbits_w_transient, SNR);

[rr, rr_sampled, gm] = matched_filter(rc, t0);

[qc_b, qc_a, qc_length] = transmitter_tf();
qc = impz(qc_b, qc_a, qc_length);
[c, b] = build_dfe_filters(qc, gm, t0, sigma2_a, N0, D, M1, M2);
[dec_sym_dfe, y] = dfe_filtering(c, b, rr_sampled, D);

rr_syms_only = y(t0_sampled+D+1:length(y));
[dec_bits, dec_syms] = QPSKdemodulator(rr_syms_only);

a_notail = a(1:Nbits/2);
sym_err_count = sum(dec_syms ~= a_notail);
Pe = sym_err_count / length(a_notail);

bits_notail = bits(1:Nbits); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);
end
