function [Pbit, Pe, err_count, sym_err_count] = simulate_viterbi_alternative(Nbits, SNR)
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
states = combvec(alphabet.',alphabet.',alphabet.',alphabet.').';    %%%%must be M2+1 times "",alphabet.'""%%%%%%%%%%

psi = psi(t0_sampled+D+1:t0_sampled+D+M2+1);
psi = repmat(psi,1,4^(M2+1));

uk = states.*psi.';

rho_k = zeros(4^(M2+1),M2+1);

distance = rho_k(:,1);
cost = zeros(4^(M2+1),length(y)-trellis_depth);
cost(:,1) = ones(4^(M2+1),1);       %initial state cost
for i = 1:length(y)-trellis_depth
    rho_k = repmat(y(i:M2+i),1,4^(M2+1)).';
    distance = (abs(rho_k-uk).^2)*ones(M2+1,1);
    [cost(:,i+1), I] = min(cost(:,i) + distance);
    if i == 1
        survivors = states(I,:);
    else
        survivors(M2+i) = states(I,M2+1);
    end
end


% Check errors
syms_tot_cut_first = survivors(t0_sampled+D+1:length(survivors)); % Cut the samples before the first symbol
sym_err_count = sum(syms_tot_cut_first.' ~= a(1:length(syms_tot_cut_first)));
Pe = sym_err_count / length(syms_tot_cut_first);

% Bit errors
dec_bits = QPSKdemodulator(syms_tot_cut_first);
bits_notail = bits(1:2*length(syms_tot_cut_first)); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);


end
