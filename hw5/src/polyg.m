function g = polyg
% input(a,b,c,d,f1,f2)
syms t
syms f1 f2 a b c d
f = t*((a*t+b)^2+f2^2*(c*t+d)^2)^2-(a*d-b*c)*(1+f1^2*t^2)^2*(a*t+b)*(c*t+d);
g = matlabFunction(flip(coeffs(f,t)));
end