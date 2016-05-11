close all; clear all; clc;

b = [zeros(1,10), 0.7424];   %beta
a = [1, -0.67];     %1-alfa

k = 0:999;
x = sin(2*pi*0.01*k);
y = filter(b,a,x);

[y1, zf] = filter(b,a,x(1:500));
[y2, zf] = filter(b,a,x(501:1000), zf);



filt = dfilt.df2(b,a);
filt.PersistentMemory = true;
y11 = filt.filter(x(1:500));
y22 = filt.filter(x(501:1000));

y_dfilt = [y11, y22];

figure;
plot(k,y);
hold on;
plot(k, [y1, y2]);
plot(k, y_dfilt);

