function detected_syms = dfe_filtering(c, b, x, D)
C = dfilt.df2(c, 1);
%C.PersistentMemory = true;
B = dfilt.df2(b, 1);
B.PersistentMemory = true;

Nsym = length(x);
detected_syms = zeros(Nsym, 1);
xff = C.filter(x);
%xfb = zeros(Nsym, 1);
%y = zeros(Nsym, 1);
for i=0:Nsym-1
     if i-D >= 0
         
        xfb_i = B.filter(detected_syms(i-D+1));
    else
        xfb_i = B.filter(0);
    end
    y_i = xff(i+1) + xfb_i;
    [~, detected_syms(i+1)] = QPSKdemodulator(y_i);
end
end
