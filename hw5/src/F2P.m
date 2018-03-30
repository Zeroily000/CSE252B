function P = F2P(wu,wv,s)
U = AAdeparameterization(wu);
V = AAdeparameterization(wv);
sigma = Xdeparameterization(s);
Z = [0 -1 0
     1  0 0
     0  0 1];
d = [sigma(1);sigma(2);(sigma(1)+sigma(2))/2];
m = U*Z*diag(d)*V.';
P = [m,-U(:,3)];
end