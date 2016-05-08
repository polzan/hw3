close all; clear all; clc;

t_0 = 33; % Must be >= qc_length (causal matched filter) ?
D = 9;

% random bits
Nbits = 1e4;
bits = round(rand(Nbits, 1));

[r_c, s_c, w, sigma2_a, N0] = transmit_bits(bits, 11, t_0 + D*4); % add padding or tx more bits?

%match filter
qc_length = 34; % after < 5e-5
[qc_b, qc_a] = transmitter_tf();
q_c = impz(qc_b, qc_a, qc_length);
q_match = flip(conj(q_c));

q_R = conv(q_c, q_match);
figure;
subplot(1,2,1);
stem(-length(q_c)+1:0, q_match);
subplot(1,2,2);
hold on;
stem(-length(q_c)+1:length(q_c)-1, q_R);
plot([t_0, t_0] - length(q_c)+1, ylim);

r_R = filter(q_match, 1, r_c); % q_match as FIR 
r_R_t0 = r_R(t_0+1:length(r_R));
r_sampled = downsample(r_R_t0, 4);

%%%c estimation%%%
r_qc = downsample(q_R, 4, mod(t_0,4));

R_QC = fft(r_qc);
C = 1./R_QC;
c = ifft(sigma2_a./(N0 + sigma2_a.*R_QC));



figure;
plot(abs(C));
figure;
stem(c);

psi = conv(r_qc, c);
figure;
stem(psi);
hold on;
plot([D+1 D+1], ylim);

received = filter(c, 1, r_sampled);
received = received(D+1:length(received));

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
bit_est = QPSKdemodulator(received);    %with c
bit_est2 = QPSKdemodulator(r_sampled(1:length(r_sampled)-D));      %without c

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
