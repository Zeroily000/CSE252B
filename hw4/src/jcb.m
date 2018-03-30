function [A,B1,B2] = jcb(h_hat, x_hat, f_A, f_B1, f_B2)
n = size(x_hat,2);
A = zeros(2,8,n);
B1 = zeros(2,2,n);
B2 = zeros(2,2,n);
h1 = h_hat(1);
h2 = h_hat(2);
h3 = h_hat(3);
h4 = h_hat(4);
h5 = h_hat(5);
h6 = h_hat(6);
h7 = h_hat(7);
h8 = h_hat(8);
for i = 1:n
    x = x_hat(1,i);
    y = x_hat(2,i);
    A(:,:,i) = double(f_A(h1, h2, h3, h4, h5, h6, h7, h8, x, y));
    B1(:,:,i) = double(f_B1(x, y));
    B2(:,:,i) = double(f_B2(h1, h2, h3, h4, h5, h6, h7, h8, x, y));
end
end