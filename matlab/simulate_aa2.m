function [Pbit, Pe, err_count, sym_err_count] = simulate_aa2(Nbits, SNR)
t0 = 20;
M1 = 3;
M2 = 2;
D = 1;
t0_sampled = floor(t0/2);
Nbits_w_transient = Nbits + t0_sampled + 2*D;
[rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(Nbits_w_transient, SNR, 1e5, 1e8);
[qc_b, qc_a, qc_length] = transmitter_tf();
qc = impz(qc_b, qc_a, qc_length);
gaa = fir1(20,0.7,'low');
rr = filter(gaa, 1, rc);
rr_sampled = downsample(rr, 2, mod(t0,2));

[c, b] = build_dfe_filters_T2(qc, gaa, t0, sigma2_a, N0, D, M1, M2);

[dec_sym_dfe, y] = dfe_filtering_2(c, b, rr_sampled, D);

rr_syms_only = y(t0_sampled/2+D+1:length(y));
[dec_bits, dec_syms] = QPSKdemodulator(rr_syms_only);

a_notail = a(1:Nbits/2);
sym_err_count = sum(dec_syms ~= a_notail);
Pe = sym_err_count / length(a_notail);

bits_notail = bits(1:length(dec_bits)); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);

end

