function x_hat0 = estimate(p_hat0,X_bar_norm)
p_bar0 = deparameterization(p_hat0);
P_hat0 = reshape(p_bar0,4,3)';

x_bar0 = P_hat0*X_bar_norm;
x_hat0 = x_bar0(1:2,:)./x_bar0(3,:);
x_hat0 = reshape(x_hat0,[],1);
end