function [bits, detected_syms] = dfe_filtering(c, b, x, D)
C = dfilt.df2(c, 1);
C.PersistentMemory = true;
B = dfilt.df2(b, 1);
B.PersistentMemory = true;

K = length(x);
detected_syms = zeros(K, 1);
xff = zeros(K, 1);
xfb = zeros(K, 1);
y = zeros(K, 1);
for i=1:K
    xff(i) = C.filter(x(i));
    if i-D >= 1
        xfb(i) = B.filter(detected_syms(i-D));
    else
        xfb(i) = B.filter(0);
    end
    y(i) = xff(i) + xfb(i);
    [bits(2*i+1:2*i+2), detected_syms(i)] = QPSKdemodulator(y(i));
end
end
