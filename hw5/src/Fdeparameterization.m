function F = Fdeparameterization(wu,wv,s)
U = AAdeparameterization(wu);
V = AAdeparameterization(wv);
sigma = Xdeparameterization(s);
Sigma = diag([sigma;0]);
F = U*Sigma*V';
end