function s = all_states(alphabet, L)
if L > 1
    s_ = all_states(alphabet, L-1);
    s = [];
    for i=1:length(alphabet)
        for j=1:size(s_,1)
            r = [s_(j,:), alphabet(i)];
            s = [s; r];
        end
    end
elseif L == 1
    for i=1:length(alphabet)
        s(i,1) = alphabet(i);
    end
else
    error('L<1');
end
end
