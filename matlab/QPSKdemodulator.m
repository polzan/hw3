function [bits, detected_syms] = QPSKdemodulator(data)
bits = zeros(2*length(data),1);
detected_syms = zeros(length(data), 1);
sym_d = 1;
v = 0;
for i = 1:length(data)
    if real(data(i)) > 0 && imag(data(i)) > 0
        bits(v+i) = 0;
        bits(v+i+1) = 0;
        detected_syms(i) = sym_d * (1 + 1j);
    elseif  real(data(i)) > 0 && imag(data(i)) < 0
        bits(v+i) = 1;
        bits(v+i+1) = 0;
        detected_syms(i) = sym_d * (1 - 1j);
    elseif real(data(i)) < 0 && imag(data(i)) > 0
        bits(v+i) = 0;
        bits(v+i+1) = 1;
        detected_syms(i) = sym_d * (-1 + 1j);
    else
        bits(v+i) = 1;
        bits(v+i+1) = 1;
        detected_syms(i) = sym_d * (-1 - 1j);
    end
    v = v+1;
end
end
