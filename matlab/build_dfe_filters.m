function [c,b] = build_dfe_filters(sigma2_a, N0, q_match_impz, h, D, M1, M2)
%assert(M1 <= D);
N1 = abs(min(h.getIndices()));
N2 = abs(max(h.getIndices()));

h_ = h.flip();
h_.setAutoPadding(true);
P = sigma2_a .* conj(h_.getData((-M1+1:0)+D));

r_w = ArrayWithIndices(N0 .* xcorr(q_match_impz), -length(q_match_impz)+1);
R = zeros(M1, M1);
h.setAutoPadding(true);
for p=0:M1-1;
    for q=0:M1-1
        s1 = (h.getData((-N1:N2) - p + q))' * h.getData(-N1:N2);
        if M2 ~= 0
            s2 = (h.getData((1:M2) +D -p))' * h.getData((1:M2) + D -q);
        else
            s2 = 0;
        end
        
        R(p+1,q+1) = sigma2_a * (s1 - s2) + r_w.getData(p-q);
    end
end
h.setAutoPadding(false);

c = linsolve(R, P);


if M2 == 0
    b = [];
    return;
end
h_D = h.translate(D);
b = zeros(M2, 1);
for i=1:M2
    b(i) = - c.' * h_D.getData(i:i+M1-1);
end
end
