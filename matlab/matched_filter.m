function [rr, rr_sampled, gm] = matched_filter(rc, t0)
[qc_b, qc_a, qc_length] = transmitter_tf();
qc = impz(qc_b, qc_a, qc_length);
gm = flip(conj(qc));
if length(gm)-1 > t0
    error('t0 is too small, gm is not causal');
end
rr = filter(gm, 1, rc);
% drop the final transient: the last bits are be lost, send more bits to
% estimate Pbit
rr_sampled = downsample(rr, 4, mod(t0,4));
% The downsampled signal's initial transient stops at floor(t0/4)
end
