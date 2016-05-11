
[rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(100, 11, 8);

figure;
hold on;
scatter(real(a), imag(a), 'x');
scatter(real(sc), imag(sc), '.');
scatter(real(rc), imag(sc), 'O');

[dec_bits, dec_syms] = QPSKdemodulator(downsample(rc, 4));

err_count = sum(abs(dec_bits-bits));
Pbit = err_count / 100;
fprintf('Pbit %f\n', Pbit);

