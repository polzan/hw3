close all; clear all; clc;

% TRANSMITTER

%freq response of Q_c(f)
T = 1;
up_factor = 4;
T_Q = T/up_factor;      % F_4 = 4

Q_c_num = [zeros(1,10), 0.7424];   %beta
Q_c_den = [1, -0.67];     %1-alfa

q_c = impz(Q_c_num, Q_c_den, 27);

[Q_f_plot, f_plot] = freqz(Q_c_num, Q_c_den, 1000, 1/(2*T_Q));

%modulator
Nbits = 1e4;
assert(mod(Nbits, 2) == 0); % Must be even
bits = round(rand(Nbits, 1));
QPSK_mod = comm.QPSKModulator('BitInput',true);     %bit grey coded
a = step(QPSK_mod,bits);
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
M1 = 19;        %at least 12
q_match = fliplr(q_c(1:M1)');  %already shifted
r_r = conv(r_c,q_match);
to = 19;  %%%%%%dal grafico di r_qc??? = E_q??(8.24) real t0= M1??
r_sampled = zeros(1,length(a_Q)/4);
for i = 1:length(a_Q)/4      %%sampling in to+i*T
    r_sampled(1,i) = r_r(to+(i-1)*up_factor*T);
end

%%%c estimation%%%
r_qc = conv(q_c,q_match);      %%real t=14??

r_qc = decimate(circshift(r_qc,length(r_qc)-to+1),4);
r_qc = circshift(r_qc,length(r_qc)/2);
R_QC = fft(r_qc);
C = 1./R_QC;
c = ifft(sigma2_a./(N0 + sigma2_a.*R_QC));
received = conv(r_sampled,c);

%decoded bits
bit_est = QPSKdemodulator(received(1,length(c):length(received)));    %with c
bit_est2 = QPSKdemodulator(r_sampled);      %without c

P_bit = (length(find(bit_est - transpose(bits))) ~= 0)/length(bits)
p_bit_no_c = length(find(bit_est2 - transpose(bits)) ~= 0)/length(bits)

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
stem(bit_est-bits');
subplot(1,2,2)
stem(bit_est2-bits');

%interf before anf after c
figure;
subplot(1,2,1)
stem(0:length(r_qc)-1,r_qc);
subplot(1,2,2)
stem(conv(c,r_qc))