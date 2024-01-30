function [w] = wn(t, x, n, T, D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
w = 0;
for k = 0:(n-1)
    for j = 0:(n-1)
        w = w + 8*sqrt(T*D) .* sin((2*k + 1) .* pi .* t ./ (2*T)) ...
            .* sin((2.*j + 1) .* pi .* x / (2.*D)) .* normrnd(0,1, size(x)) ... 
            ./ (pi^2 .* (2*k+1).*(2.*j + 1));
    end
end
end