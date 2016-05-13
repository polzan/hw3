function [b, a, qc_length, qc] = transmitter_tf()
b = [zeros(1,10), 0.7424];   %beta
a = [1, -0.67];     %1-alfa
qc_length = 34; % After this < 5e-5
qc = impz(b, a, qc_length);
end
