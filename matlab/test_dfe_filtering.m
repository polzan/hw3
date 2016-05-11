close all; clear all; clc;

c = [1];
b = [-0.5];

K = 1e3;
bits = round(rand(2*K, 1));
x = QPSKmodulator(bits);
D = 2;

[bits, a] = dfe_filtering(c,b,x,D);

figure;
hold on;
stem(100:119, real(x(101:120)));
stem((100:119), real(a(101:120)));
