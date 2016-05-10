function [b,a, qc_length] = transmitter_tf()
b = [zeros(1,10), 0.7424];   %beta
a = [1, -0.67];     %1-alfa
qc_length = 34;
end
