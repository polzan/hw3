function [r_c, s_c, w] = transmit_bits(bits, SNR, padding, varargin)
zeronoise = false;
if ~isempty(varargin)
    if any(strcmp(varargin, 'zeronoise'))
        zeronoise = true;
    end
end

T = 1;
up_factor = 4;
T_Q = T/up_factor;      % F_Q = 4/T

qc_length = 32;
[qc_b, qc_a] = transmitter_tf();
q_c = impz(qc_b, qc_a, qc_length);

% bitmap
Nbits = length(bits);
assert(mod(Nbits, 2) == 0); % Must be even
a = QPSKmodulator(bits);
a_Q = upsample(a,up_factor);

%compute additional params.
E_q = norm(q_c)^2;      %energy of q
sigma2_a = 1/sqrt(2);%/up_factor;    %%%%%%%%%%%????????
%%%%lesson 21/04??
sigma2_w = 10*log10(E_q*sigma2_a) - SNR;       %sigma2 = Msc*Q0 - SNR in dB
sigma2_w = 10^(sigma2_w/10);            %%%linear noise power
N0 = sigma2_w*T_Q;

% noise
w = sqrt(sigma2_w).*(randn(length(a_Q) + padding, 1) + 1j*(randn(length(a_Q) + padding, 1)));

%received signal
s_c = filter(qc_b, qc_a, [a_Q; zeros(padding, 1)]);
if zeronoise
    r_c = s_c;
else
    r_c = s_c + w;
end
end
