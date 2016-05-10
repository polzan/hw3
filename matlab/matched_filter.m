function q_match = matched_filter()
[qc_b, qc_a, qc_length] = transmitter_tf();
q_c = impz(qc_b, qc_a, qc_length);
q_match = flip(conj(q_c));
end
