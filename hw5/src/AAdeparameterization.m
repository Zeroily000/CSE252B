function UV = AAdeparameterization(w)
theta = norm(w);
wx = skew(w);
UV = cos(theta)*eye(3) + sinc(theta/pi)*wx + (1-cos(theta))*(w*w.')/theta^2;
end