close all; clear all; clc;

% TRANSMITTER

%freq response of Q_c(f)
T = 1;
up_factor = 4;
T_Q = T/up_factor;      % F_4 = 4

Q_c_num = [zeros(1,10), 0.7424];   %beta
Q_c_den = [1, -0.67];     %1-alfa

q_c = impz(Q_c_num, Q_c_den, 32);

[Q_f_plot, f_plot] = freqz(Q_c_num, Q_c_den, 1000, 1/(2*T_Q));

%modulator
Nbits = 1e4;
assert(mod(Nbits, 2) == 0); % Must be even
bits = round(rand(Nbits, 1));
a = QPSKmodulator(bits);
a_Q = upsample(a,up_factor);

%compute additional params.
SNR_target = 11;
E_q = norm(q_c)^2;      %energy of q
sigma2_a = 1/sqrt(2);%/up_factor;    %%%%%%%%%%%????????
%%%%lesson 21/04??
sigma2_w = 10*log10(E_q*sigma2_a) - SNR_target;       %sigma2 = Msc*Q0 - SNR in dB
sigma2_w = 10^(sigma2_w/10);            %%%linear noise power
N0 = sigma2_w*T_Q;
w = sqrt(sigma2_w).*(randn(length(a_Q),1) + 1j*(randn(length(a_Q),1)));

%rc signal
s_c = filter(Q_c_num, Q_c_den, a_Q);
r_c = s_c + w;

% RECEIVER

%match filter
q_match = flip(conj(q_c));
%[~, t_0_i] = max(abs(q_match));
%t_0 = t_0_i-1;
t_0 = 32;

% figure;
% Q_match = fft(q_match, 2048);
% Q_f = fft(q_c, 2048);
% f = linspace(0, 1/T_Q, 2048);
% hold on;
% plot(f, imag(Q_f));
% plot(f, imag(Q_match));

q_R = conv(q_c, q_match);

figure;
hold on;
stem(q_R);
plot([t_0, t_0], ylim);

r_R = filter(q_match, 1, r_c); % q_match as FIR 

r_R_t0 = r_R(t_0+1:length(r_R));

r_sampled = downsample(r_R_t0, 4);

% figure;
% hold on;
% stem(0:length(r_R)-1, abs(r_R));
% stem(t_0:length(r_R)-1, abs(r_R_t0));
% stem(t_0:4:length(r_R)-1, abs(r_sampled));

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
fprintf('Pbit = %f, %f (no C)\n', P_bit, p_bit_no_c);

figure;
subplot(1,3,1);
plot(f_plot, 20*log10(abs(Q_f_plot)));    %2/T but qc works T/4 >> plot half a period
xlabel('n/T_Q'); ylabel('|Q_c(f)|');
subplot(1,3,2);
stem(0:length(q_c)-1, q_c);
xlabel('nT'); ylabel('q_c(nT)');
subplot(1,3,3)
stem(0:length(q_match)-1, q_match);

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