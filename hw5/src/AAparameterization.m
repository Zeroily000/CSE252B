function w = AAparameterization(UV)
[~,~,V] = svd(UV-eye(3));
v = V(:,end);
v_hat = zeros(3,1);
v_hat(1) = UV(3,2) - UV(2,3);
v_hat(2) = UV(1,3) - UV(3,1);
v_hat(3) = UV(2,1) - UV(1,2);
sin_theta = v.'*v_hat/2;
cos_theta = (trace(UV) - 1)/2;
theta = atan2(sin_theta,cos_theta);
w = theta * v / norm(v);
w = w*(1-2*pi/theta*ceil((theta-pi)/(2*pi)));
end