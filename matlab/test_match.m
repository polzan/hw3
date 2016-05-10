close all;
[rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(100, 11, 1e5, 1e8);

t0 = 33;

[rr, rr_sampled, gm] = matched_filter(rc, t0);
[qc_b, qc_a, qc_length] = transmitter_tf();
qc = impz(qc_b, qc_a, qc_length);

figure;
subplot(1,2,1);
stem((-length(gm)+1:0) + t0, gm);
subplot(1,2,2);
qR = conv(qc, gm);
stem((-length(qR)+1:0) + t0, qR);

figure;
hold on;
scatter(real(a), imag(a), 'x');
%scatter(real(sc), imag(sc), '.');
%scatter(real(rc), imag(sc), 'O');
scatter(real(rr), imag(rr), '.');

figure;
subplot(2,1,1);
hold on;
stem((0:length(rr)-1), real(rr));
plot(t0*ones(2,1), ylim);
subplot(2,1,2);
hold on;
dropped = mod(t0,4);
t0_sampled = floor(t0/4);
stem((0:length(rr_sampled)-1), real(rr_sampled));
plot(t0_sampled*ones(2,1), ylim);

rr_syms_only = rr_sampled(t0_sampled+1:length(rr_sampled));
[dec_bits, dec_syms] = QPSKdemodulator(rr_syms_only);

figure;
hold on;
scatter(real(rr_syms_only), imag(rr_syms_only));
scatter(real(a), imag(a), 'x');

bits_notail = bits(1:length(dec_bits)); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);
fprintf('Pbit %f\n', Pbit);

