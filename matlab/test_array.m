close all; clear all; clc;

[b,a] = transmitter_tf();

h = ArrayWithIndices(impz(b,a), -10);

%h.alignOffset(-1);

g = h.flip();

figure;
hold on;
stem(h.getIndices(), h.getAll());
figure;
h.setAutoPadding(true);
stem(-50:50, h.getData(-50:50));

figure;
stem(g.getIndices(), g.getAll());

h.setData(-50, 1)
h.setData(-40:-36, -1*ones( 5, 1));
figure;
hold on;
stem(h.getIndices(), h.getAll());

h = ArrayWithIndices(impz(b,a), 0);
g = h.flip();
% g = ArrayWithIndices();
% g.setData(0, 1);
% g.alignOffset(-100);

c = conv_indices(h, g);

figure;
stem(c.getIndices(), c.getAll());




