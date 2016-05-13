function a = QPSKmodulator(bits)
if mod(length(bits), 2) ~= 0
    error('The QPSK modulator needs an even number of bits');
end
if size(bits, 1) == 1
    bits = transpose(bits);
end
a = zeros(length(bits)/2, 1);
d = 1;
for i=2:2:length(bits)
    if bits(i)
        a(i/2) = -d;
    else
        a(i/2) = d;
    end
    if bits(i-1)
        a(i/2) = a(i/2) - 1j*d;
    else
        a(i/2) = a(i/2) + 1j*d;
    end
end
end
