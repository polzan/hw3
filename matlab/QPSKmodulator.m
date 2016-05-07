function a = QPSKmodulator(bits)
if size(bits, 1) == 1
    bits = transpose(bits);
end
QPSK_mod = comm.QPSKModulator('BitInput',true);     %bit grey coded
a = step(QPSK_mod,bits);
end
