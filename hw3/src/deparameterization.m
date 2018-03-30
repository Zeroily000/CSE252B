function R = deparameterization(w)
theta = norm(w);
wx = skew(w);
R = cos(theta)*eye(3) + sinc(theta/pi)*wx + (1-cos(theta))/theta^2*w*w';
end