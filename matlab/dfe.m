close all;

t_0 = 33; % Must be >= qc_length (causal matched filter) ?
D = 9;

% random bits
Nbits = 1e4;
bits = round(rand(Nbits, 1));

bits(1:500) = ones(500, 1);
bits(501:1000) = zeros(500, 1);

SNR_target = 11;   %lin. scale
[r_c, s_c, w, sigma2_a, N0] = transmit_bits(bits, SNR_target, t_0 + D*4, 'zeronoise'); % add padding or tx more bits?

%match filter
q_match = matched_filter();

q_R = conv(q_c, q_match);
figure;
subplot(1,2,1);
stem(-length(q_c)+1:0, q_match);
subplot(1,2,2);
hold on;
stem(-length(q_c)+1:length(q_c)-1, q_R);
plot([t_0, t_0] - length(q_c)+1, ylim);

r_R = filter(q_match, 1, r_c); % q_match as FIR 
r_sampled = downsample(r_R, 4, mod(t_0, 4));

%%%c estimation%%%
r_qc = downsample(q_R, 4, mod(t_0,4));
t_0_sampled = floor(t_0 / 4);

r_qc_t0 = ArrayWithIndices(r_qc, -t_0_sampled);


figure;
stem(r_qc_t0.getIndices(), r_qc_t0.getAll());

M1 = 7;
N2 = length(r_qc) - t_0_sampled - 1;
M2 = -D + N2 + M1 -1;
[c, b] = dfe_filters(sigma2_a, N0, downsample(q_match, 4, mod(t_0, 4)), r_qc_t0, D, M1, M2);

figure;
stem(c);

psi = conv(r_qc, c);
figure;
stem(psi);
hold on;
plot([D+1 D+1], ylim);

received = dfe_filtering(c,b,r_sampled,D);


skip_transient = t_0_sampled + D;
received = received(skip_transient+1:length(received));

% c_opt_b = [zeros(D,1); sigma2_a];
% c_opt_a = sigma2_a .* r_qc;
% c_opt_a(1) = c_opt_a(1) + N0;



%figure;
%stem(r_qc);

%figure;
%impz(c_opt_b, c_opt_a);

%psi_b = conv(r_qc, c_opt_b);
%psi_a = c_opt_a;

%figure;
%impz(psi_b, psi_a);

%received = filter(c_opt_b, c_opt_a, r_sampled);

%decoded bits
bit_est = QPSKdemodulator([received(3:length(received)); 0; 1; 0]);    %with c
bit_est2 = QPSKdemodulator(r_sampled(t_0_sampled + 1 :length(r_sampled)-D));      %without c

diff_bits1 = bit_est - bits(1:length(bit_est));
err_count1 = sum(abs(diff_bits1 ~=0));

diff_bits2 = bit_est2 - bits(1:length(bit_est2));
err_count2 = sum(abs(diff_bits2 ~=0));

P_bit = err_count1 / length(bit_est);
p_bit_no_c = err_count2 / length(bit_est2);
fprintf('Pbit = %f, %f (no C)\n', P_bit, p_bit_no_c);

%compares error with and without c (something wrong)
figure;
subplot(1,2,1)
stem(bit_est-bits);
subplot(1,2,2)
stem(bit_est2-bits);

%interf before anf after c
figure;
subplot(1,2,1)
stem(0:length(r_qc)-1,r_qc);
subplot(1,2,2) 
stem(conv(c,r_qc))


figure;
subplot(3,1,1);
stem(0:length(bits)-1, bits);
subplot(3,1,2);
stem(t_0_sampled+(0:length(bit_est2)-1), bit_est2);
subplot(3,1,3);
stem((0:length(bit_est)-1) + t_0_sampled + D, bit_est);