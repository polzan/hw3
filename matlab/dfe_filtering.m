function [detected_syms, y] = dfe_filtering(c, b, x, D)
xff = filter(c, 1, x);
if isempty(b) || all(b == 0)
    y = xff;
    [~, detected_syms] = QPSKdemodulator(y);
else
    Nsym = length(x);
    detected_syms = zeros(Nsym, 1);
    y = zeros(Nsym, 1);
    for k=0:Nsym-1
        if k-D-length(b) >= 0
            a_D = detected_syms((k-D-length(b)+1:k-D)+1);
        elseif k-D >= 0
            a_D = [zeros(length(b)-k+D-1, 1); detected_syms((0:k-D)+1)];
        else
            a_D = zeros(length(b), 1);
        end
        xfb_k = transpose(b) * flip(a_D);
        y(k+1) = xff(k+1) + xfb_k;
        [~, detected_syms(k+1)] = QPSKdemodulator(y(k+1));
    end
end
end
