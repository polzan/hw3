function c = conv_indices(h, g)
minoff = min(h.getOffset(), g.getOffset());
h.alignOffset(minoff);
g.alignOffset(minoff);
c_data = conv(h.getAll(), g.getAll());
c_off = 2*minoff;
c = ArrayWithIndices(c_data, c_off);
end

