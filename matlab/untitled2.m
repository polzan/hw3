close all; clear all; clc;

%freq response of Q_c(f)
T = 1;
up_factor = 4;
T_4 = T/up_factor;      % F_4 = 4
b = zeros(1,11);
b(1,11) = 0.7424;   %beta
a = [1, -0.67];     %1-alfa
Q_f = freqz(b,a,'whole',2000);
q = ifft(Q_f);
q = q(1:27);            %works at T_4
f = linspace(0,2/T,1000);


%modulator
D = 13;
L = 4095;
PN = rand(1,3000)>0.5;  %pn_seq(L,1,D);        %must be even length for the following line
QPSK_mod = comm.QPSKModulator('BitInput',true);     %bit grey coded
a_k = step(QPSK_mod,PN');
a_k = upsample(a_k,up_factor);


%compute additional params.
SNR_target = 11;
E_q = q.'*conj(q);      %energy of q
sigma2_a = 1/sqrt(2);%/up_factor;    %%%%%%%%%%%????????
%%%%lesson 21/04??
sigma2_omega = 10*log10(E_q*sigma2_a) - SNR_target;       %sigma2 = Msc*Q0 - SNR in dB
sigma2_omega = 10^(sigma2_omega/10);            %%%linear noise power
N0 = sigma2_omega*T_4;
omega = sqrt(sigma2_omega).*(randn(length(a_k)+length(q)-1,1) + 1j*(randn(length(a_k)+length(q)-1,1)));

%rc signal
r_c = conv(a_k,q) + omega;


%match filter
M1 = 19;        %at least 12
q_match = fliplr(q(1:M1)');  %already shifted
r_r = conv(r_c,q_match);
to = 19;  %%%%%%dal grafico di r_qc??? = E_q??(8.24) real t0= M1??
r_sampled = zeros(1,length(a_k)/4);
for i = 1:length(a_k)/4      %%sampling in to+i*T
    r_sampled(1,i) = r_r(to+(i-1)*up_factor*T);
end

%%%c estimation%%%
r_qc = conv(q,q_match);      %%real t=14??

r_qc = decimate(circshift(r_qc,length(r_qc)-to+1),4);
r_qc = circshift(r_qc,length(r_qc)/2);
R_QC = fft(r_qc);
C = 1./R_QC;
c = ifft(sigma2_a./(N0 + sigma2_a.*R_QC));
received = conv(r_sampled,c);

%decoded bits
bit_est = QPSKdemodulator(received(1,length(c):length(received)));    %with c
bit_est2 = QPSKdemodulator(r_sampled);      %without c

P_bit = length(find(bit_est - PN ~= 0))/length(PN)
p_bit_no_c = length(find(bit_est2 -PN ~= 0))/length(PN)

figure;
subplot(1,3,1);
plot(f, 10*log10(abs(Q_f(1:1000))));    %2/T but qc works T/4 >> plot half a period
xlabel('n/T'); ylabel('|Q_c(f)|');
subplot(1,3,2);
stem(0:length(q)-1, q);
xlabel('nT'); ylabel('q_c(nT)');
subplot(1,3,3)
stem(0:length(q_match)-1, q_match);

%compares error with and without c (something wrong)
figure(2);
subplot(1,2,1)
stem(bit_est-PN);
subplot(1,2,2)
stem(bit_est2-PN);

%interf before anf after c
figure(3);
subplot(1,2,1)
stem(0:length(r_qc)-1,r_qc);
subplot(1,2,2)
stem(conv(c,r_qc))