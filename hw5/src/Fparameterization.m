function [wu,wv,s] = Fparameterization(F)
[U,Sigma,V] = svd(F);
if round(det(U)) == -1
    U = -U;
end
if round(det(V)) == -1
    V = -V;
end
sigma = [Sigma(1,1);Sigma(2,2)];
wu = AAparameterization(U);
wv = AAparameterization(V);
s = Xparameterization(sigma);
end