function [Pbit, Pe, err_count, sym_err_count] = simulate_viterbi(Nbits, SNR)
% Parameters
t0 = 33;
t0_sampled = floor(t0/4);
M1 = 3;
M2 = 2;
D = 1;

% Viterbi params
trellis_depth = 20;
L1 = 0; % No precursors
L2 = M2; % Same postcursors of the DFE

Nbits_w_transient = Nbits + 2*t0_sampled + 2*D + 2*(trellis_depth-1);

% Generate a transmitted signal with random iid bits
[rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(Nbits_w_transient, SNR);

% Matched filter + sampling at t0
[rr, rr_sampled, gm] = matched_filter(rc, t0);

% Use only the feedforward filter of a DFE
[~, ~, ~, qc] = transmitter_tf();
[c, ~] = build_dfe_filters(qc, gm, t0, sigma2_a, N0, D, M1, M2);
y = filter(c, 1, rr_sampled);

% Compute the overall imp.resp. psi (w. total delay t0_sampled + D)
qR = conv(qc, gm);
qR_sampled = downsample(qR, 4, mod(t0, 4));
psi = conv(qR_sampled, c);

alphabet = [1+1j; 1-1j; -1-1j; -1+1j];
vd = ViterbiDetector(alphabet, trellis_depth, L1, L2, psi, t0_sampled + D);
detected_syms = vd.detect(y);

% Check errors
syms_tot_cut_first = detected_syms(t0_sampled+D+(trellis_depth-1)+1:length(detected_syms)); % Cut the samples before the first symbol
sym_err_count = sum(syms_tot_cut_first ~= a(1:(Nbits/2)));
Pe = sym_err_count / (Nbits/2);

% Bit errors
dec_bits = QPSKdemodulator(syms_tot_cut_first);
bits_notail = bits(1:Nbits); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);
end
