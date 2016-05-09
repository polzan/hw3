close all; clear all; clc;

[b,a] = transmitter_tf();

h = ArrayWithIndices(impz(b,a));

h.offset = -5;

g = h.flip();

figure;
hold on;
stem(h.indices, h.all);
h.autoPadding = true;
stem(-50:50, h(-50:50));

figure;
stem(g.indices, g.all);
