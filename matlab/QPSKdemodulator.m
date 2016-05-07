function [ a ] = QPSKdemodulator(data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
a = zeros(2*length(data),1);
v = 0;
for i = 1:length(data)
    if real(data(i)) > 0 && imag(data(i)) > 0
        a(v+i) = 0;
        a(v+i+1) = 0;
    elseif  real(data(i)) > 0 && imag(data(i)) < 0
        a(v+i) = 1;
        a(v+i+1) = 0;
    elseif real(data(i)) < 0 && imag(data(i)) > 0
        a(v+i) = 0;
        a(v+i+1) = 1;
    else
       a(v+i) = 1;
        a(v+i+1) = 1;
    end
    v = v+1;
end

end

