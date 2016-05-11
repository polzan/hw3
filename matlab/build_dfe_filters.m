function [c,b] = build_dfe_filters(qc, gm, t0, sigma2_a, N0, D, M1, M2)
h_T_Q = conv(qc, gm);
h = downsample(h_T_Q, 4, mod(t0, 4));
t0_sampled = floor(t0/4);
N1 = t0_sampled;
N2 = length(h)-1 - t0_sampled;
assert(length(h) == N1 + N2 + 1);
if nargin < 8
    M2 = -D + N2 + M1 -1;
end
r_gm = xcorr(gm, 'none');
r_gm_zero_offset = (length(r_gm)-1)/2;
r_gm_sampled = downsample(r_gm, 4, mod(r_gm_zero_offset, 4));
r_gm_sampled_zero_offset = floor(r_gm_zero_offset / 4);

r_ax = sigma2_a .* conj(flip(h));
r_h = xcorr(h, 'none');
r_h_zero_offset = (length(r_h)-1)/2;
r_w = N0 .* r_gm_sampled;
r_w_zero_offset = r_gm_sampled_zero_offset;

assert(M1-D <= N1);
assert(D >= 0);
assert(D <= N2-1);
P = r_ax((1:M1) - D + t0_sampled);

%assert(D+M2 <= N2-1);
R = zeros(M1, M1);
for p=0:M1-1;
    for q=0:M1-1
        s1 = r_h(r_h_zero_offset + p-q +1);
        if M2 == 0
            s2 = 0;
        else
            pad_amount = M2;
            h_pad = [zeros(pad_amount, 1); h; zeros(pad_amount, 1)];
            h1 = h_pad((D+1:D+M2) -q + t0_sampled +1 +pad_amount);
            h2 = conj(h_pad((D+1:D+M2) -p + t0_sampled +1 +pad_amount));
            s2 = transpose(h1) * h2;
        end
        
        R(p+1,q+1) = sigma2_a * (s1 - s2) + r_w(p-q+1+r_w_zero_offset);
    end
end

c = linsolve(R, P);

if M2 == 0
    b = [0];
else
    psi = conv(h, c);
    b = zeros(M2+1, 1);
    b((1:M2)+1)= -psi((1:M2)+1+t0_sampled+D);
    b = b(2:M2+1);
end
end
