function [ I, s ] = intSimpson(f, a, b, n)
%INTSIMPSON Calculate simpson 1/3 numerical integral of function handle f
%   input function handel f and spanning a to b with an even step size n

h = (b - a)/n;
s = f(a) + f(b);

for i = 1:2:n-1
    s = s + 4*f(a+i*h);
end
for i = 2:2:n-2
    s = s + 2*f(a+i*h);
end
I = h/3 * s;
end


