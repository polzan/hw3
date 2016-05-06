function [pn] = pn_seq(L, rep, off)
    r = log2(L+1);   %order of shift register
%     seq_size = 2^r-1;       %period of the sequence
%      rep = mod(L,seq_size);
    sr = ones(1,r);     %should be different from [0 0 0 0]
    pn = zeros(1,L);
    for i = 1 : L
        if r == 8
            pn(i) = xor(sr(r-6),xor(xor(sr(r),sr(r-4)),sr(r-5)));
        elseif r == 9
            pn(i) = xor(sr(r),sr(r-4));
        elseif r == 12
            pn(i) = xor(sr(r-10),xor(xor(sr(r),sr(r-1)),sr(r-2)));
        elseif r == 13 || r == 14
            pn(i) = xor(sr(r-12),xor(xor(sr(r),sr(r-1)),sr(r-2)));
        elseif r == 16
            pn(i) = xor(sr(r-5),xor(xor(sr(r),sr(r-2)),sr(r-3)));
        elseif r == 18
            pn(i) = xor(sr(r),sr(r-7));
        elseif r == 10 || r == 17 || r == 20
            pn(i) = xor(sr(r),sr(r-3));
        elseif r == 19
            pn(i) = xor(sr(r-5),xor(xor(sr(r),sr(r-1)),sr(r-2)));
        elseif r == 5 || r == 11
            pn(i) = xor(sr(r),sr(r-2));
        else
        pn(i) = xor(sr(r),sr(r-1));  
        end
        sr(2:r) = sr(1:r-1);    %shift right
        sr(1) = pn(i);      % update first
    end
 %    pn(1,find(pn==0)) = -1;   %should add an option
%     if rep ~= 0
%     pn = [pn pn(1:1+rep-1)];  %%%%if L is not = 2^r-1 i repeat
%     end
if off ~=0
    pn = [repmat(pn,1,rep) pn(1:off)];
end
end
