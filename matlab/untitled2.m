close all; clear all; clc;

t_0 = 31;

% random bits
Nbits = 1e4;
bits = round(rand(Nbits, 1));

r_c = transmit_bits(bits, 11, t_0); % ? pad with t_0 ?

%match filter
qc_length = 32;
[qc_b, qc_a] = transmitter_tf();
q_c = impz(qc_b, qc_a, qc_length);
q_match = flip(conj(q_c));

%t_0 = 31;

q_R = conv(q_c, q_match);
figure;
hold on;
stem(0:length(q_R)-1, q_R);
plot([t_0, t_0], ylim);

r_R = filter(q_match, 1, r_c); % q_match as FIR 
r_R_t0 = r_R(t_0+1:length(r_R));
r_sampled = downsample(r_R_t0, 4);

%%%c estimation%%%
r_qc = conv(q_c,q_match);      %%real t=14??

r_qc = decimate(circshift(r_qc,length(r_qc)-t_0+1),4);
r_qc = circshift(r_qc,length(r_qc)/2);
R_QC = fft(r_qc);
C = 1./R_QC;
c = ifft(sigma2_a./(N0 + sigma2_a.*R_QC));
received = conv(r_sampled,c);

%decoded bits
bit_est = QPSKdemodulator(received(1,length(c):length(received)));    %with c
bit_est2 = QPSKdemodulator(r_sampled);      %without c

diff_bits1 = bit_est - bits(1:length(bit_est));
err_count1 = sum(abs(diff_bits1 ~=0));

diff_bits2 = bit_est2 - bits(1:length(bit_est2));
err_count2 = sum(abs(diff_bits2 ~=0));

P_bit = err_count1 / length(bit_est);
p_bit_no_c = err_count2 / length(bit_est2);
% fprintf('Pbit = %f, %f (no C)\n', P_bit, p_bit_no_c);

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
