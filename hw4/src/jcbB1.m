function f_B1 = jcbB1
syms x y
x_hat = [x y].';
x_bar = deparameterization(x_hat);
x1_hat = x_bar;
x1_hat_inhomo = x1_hat(1:2,:) ./ x1_hat(3,:);
f_B1([x y]) = jacobian(x1_hat_inhomo,x_hat);
end