function [f_A,f_B1,f_B2] = jcbFunction
syms wu1 wu2 wu3 wv1 wv2 wv3 s X_hat1 X_hat2 X_hat3
%% P1 P2
wu = [wu1;wu2;wu3];
wv = [wv1;wv2;wv3];
P2 = F2P(wu,wv,s);
P1 = [eye(3),zeros(3,1)];
%% X
X_hat = [X_hat1; X_hat2; X_hat3];
X = Xdeparameterization(X_hat);
%% x1 x2
x2 = P2*X;
x2_inhomo = x2(1:2) / x2(3);
x1 = P1*X;
x1_inhomo = x1(1:2) / x1(3);
%% jacobian
% input(X_hat1,X_hat2,X_hat3,s,wu1,wu2,wu3,wv1,wv2,wv3)
f_A = matlabFunction(jacobian(x2_inhomo,[wu.',wv.',s]));
% input(X_hat1,X_hat2,X_hat3)
f_B1 = matlabFunction(jacobian(x1_inhomo,X_hat));
% input(X_hat1,X_hat2,X_hat3,s,wu1,wu2,wu3,wv1,wv2,wv3)
f_B2 = matlabFunction(jacobian(x2_inhomo,X_hat));
end