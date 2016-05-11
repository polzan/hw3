close all;
[rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(1e4, 11, 1e5, 1e8);

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








% equalizer = dfe(M1, M2, rls(0.3));
% %equalizer.RefTap = D;
% y = equalize(equalizer, rr_sampled);
M1 = 4;
M2 = 3;
D = 4;
[c, b] = build_dfe_filters(qc, gm, t0, sigma2_a, N0, D, M1, M2);

[dec_sym_dfe, y] = dfe_filtering(c, b, rr_sampled, D);

% loop_tf_b = c;
% loop_tf_a = -b;
% loop_tf_a(1) = loop_tf_a(1) + 1;
% figure;
% impz(loop_tf_b, loop_tf_a);
% 
% y = filter(loop_tf_b, loop_tf_a, rr_sampled);

figure;
stem(c);

figure;
qR_sampled = downsample(qR, 4, mod(t0, 4));
psi = conv(c, qR_sampled);
stem(psi);

%rr_syms_only = rr_sampled(t0_sampled+1:length(rr_sampled));
rr_syms_only = y(t0_sampled+D+1:length(y));
[dec_bits, dec_syms] = QPSKdemodulator(rr_syms_only);

figure;
hold on;
%scatter(real(rr_syms_only), imag(rr_syms_only), '.');
scatter(real(rr_sampled), imag(rr_sampled), '.');
scatter(real(y), imag(y), 'O');
scatter(real(a), imag(a), 'x');

bits_notail = bits(1:length(dec_bits)); % bits transmitted without the final transient
err_count = sum(abs(dec_bits-bits_notail));
Pbit = err_count / length(bits_notail);
fprintf('Pbit %f\n', Pbit);

% figure;
% c = equalizer.Weights(1:M1);
% b = equalizer.Weights(M1+1:M2);
% hold on;
% stem(c);
% %stem(b);

x = [(1+1j) .* ones(1000, 1); (1-1j).*ones(1000,1)];
x = [x;x];

y = dfe_filtering(c,b,x,D);

figure;
hold on;
stem(imag(x));
stem(imag(y));

figure;
hold on;
scatter(real(x), imag(x));
scatter(real(y), imag(y));

loop_tf_b = c;
loop_tf_a = -b;
loop_tf_a(1) = loop_tf_a(1) + 1;
figure;
impz(loop_tf_b, loop_tf_a);

y = filter(loop_tf_b, loop_tf_a, rr_sampled);

figure;
hold on;
scatter(real(y), imag(y));
scatter(real(a), imag(a));