function [detected_syms, y] = dfe_filtering(c, b, x, D)
xff = filter(c, 1, x);
if isempty(b) || all(b == 0)
    y = xff;
    [~, detected_syms] = QPSKdemodulator(y);
else
    Nsym = length(x);
    detected_syms = zeros(Nsym, 1);
    y = zeros(Nsym, 1);
    %B = dfilt.df2(b, 1);
    %B.PersistentMemory = true;
    b_state = [];
    for k=0:Nsym-1
        if k-(D+10) >= 0
            [xfb_k, b_state] = filter(b, 1, detected_syms(k-(D+10)+1), b_state);
        else
            [xfb_k, b_state] = filter(b, 1, 0, b_state);
        end
        y(k+1) = xff(k+1) + xfb_k;
        [~, detected_syms(k+1)] = QPSKdemodulator(y(k+1));
    end
end
end
