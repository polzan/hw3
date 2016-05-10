function [rc, sc, a, wc, sigma2_a, N0] = QPSKtransmitter(bits, SNR, noise_seed)
if nargin < 3
    noise_seed = 'shuffle';
end

T = 1; % Sym period
up_factor = 4;
T_Q = T/up_factor;      % F_Q = 4/T

[qc_b, qc_a, qc_length] = transmitter_tf();
qc = impz(qc_b, qc_a, qc_length);

% bitmap
[a, abs_a] = QPSKmodulator(bits);
a_up = upsample(a,up_factor);

%compute additional params.
E_qc = norm(qc)^2;      %energy of q_c
sigma2_a = abs_a^2; % a uniform with mean 0

sigma2_wc = (E_qc * sigma2_a) / 10^(SNR/10);
N0 = sigma2_wc * T_Q;

% noise
rng_oldstate = rng(noise_seed);
wc = sqrt(sigma2_wc).*(randn(length(a_up), 1) + 1j*(randn(length(a_up), 1)));
rng(rng_oldstate);

%received signal
sc = filter(qc_b, qc_a, a_up);
rc = sc + wc;
end
