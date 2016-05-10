function [rc, sc, a, bits, wc, sigma2_a, N0] = QPSKtransmitter_random(Nbits, SNR, bits_seed, noise_seed)
if nargin < 3
    bits_seed = 'shuffle';
end
if nargin < 4
    noise_seed = 'shuffle';
end
rng_oldstate = rng(bits_seed);
bits = round(rand(Nbits, 1));
rng(rng_oldstate);

[rc, sc, a, wc, sigma2_a, N0] = QPSKtransmitter(bits, SNR, noise_seed);
end
