function [Pbit, Pe, err_count, sym_err_count] = simulate_viterbi(Nbits, SNR)
% Parameters
t0 = 33;
t0_sampled = floor(t0/4);
M1 = 3;
M2 = 3;
D = 1;

% Viterbi params
trellis_depth = 30;
L1 = 0; % No precursors
L2 = M2; % Same postcursors of the DFE

Nbits_w_transient = Nbits + 2*t0_sampled + 2*D;

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

% Run the viterbi algorithm on blocks of samples
block_num = floor(length(y)/trellis_depth);
last_block_length = mod(length(y), trellis_depth);

detected_syms = [];
initial_path_metrics = zeros(4^(L1+L2), 1); % Start with equal likelyhood for every state
for i=0:block_num-1
    rho_block = y(trellis_depth*i+(1:trellis_depth));
    [detected_syms_block, final_path_metrics] = vd.detect_symbols(rho_block, initial_path_metrics);
    detected_syms = [detected_syms; detected_syms_block];
    initial_path_metrics = final_path_metrics;
end

%Last block
if last_block_length > 0
    vd_last = ViterbiDetector(alphabet, last_block_length, L1, L2, psi, t0_sampled + D);
    rho_last = y(trellis_depth*block_num + (1:last_block_length));
    [detected_syms_last, final_path_metrics] = vd_last.detect_symbols(rho_last, initial_path_metrics);
    detected_syms = [detected_syms; detected_syms_last];
end

% Check errors
syms_tot_cut_first = detected_syms(t0_sampled+D+1:length(detected_syms)); % Cut the samples before the first symbol
sym_err_count = sum(syms_tot_cut_first ~= a(1:(Nbits/2)));
Pe = sym_err_count / (Nbits/2);

% Bit errors
dec_bits = QPSKdemodulator(syms_tot_cut_first);
bits_notail = bits(1:Nbits); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);
end
