function [detected_syms, y] = dfe_filtering(c, b, x, D)
xff = filter(c, 1, x);
if isempty(b) || all(b == 0)
    y = xff;
    [~, detected_syms] = QPSKdemodulator(y);
else
    Nsym = length(x);
    detected_syms = zeros(Nsym, 1);
    y = zeros(Nsym, 1);
    B = dfilt.df2(b, 1);
    B.PersistentMemory = true;
    for k=0:Nsym-1
        if k-D >= 0
            xfb_k = B.filter(detected_syms(k-D+1));
        else
            xfb_k = B.filter(0);
        end
        y(k+1) = xff(k+1) + xfb_k;
        [~, detected_syms(k+1)] = QPSKdemodulator(y(k+1));
    end
end
end
