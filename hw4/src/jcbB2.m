function f_B2 = jcbB2
syms h1 h2 h3 h4 h5 h6 h7 h8 x y
h_hat = [h1 h2 h3 h4 h5 h6 h7 h8].';
x_hat = [x y].';
H2 = reshape(deparameterization(h_hat),3,3).';
x_bar = deparameterization(x_hat);
x2_hat = H2 * x_bar;
x2_hat_inhomo = x2_hat(1:2,:) ./ x2_hat(3,:);
f_B2([h1 h2 h3 h4 h5 h6 h7 h8 x y]) = jacobian(x2_hat_inhomo,x_hat);
end