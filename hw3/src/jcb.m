function J = jcb(f_large,f_small,X_wrd_inhomo,w,t)
n = size(X_wrd_inhomo,2);
J = zeros(2*n,6);
theta = norm(w);
if theta < pi/36
    f = f_small;
else
    f = f_large;
end
for i = 1:n
    J(2*i-1 : 2*i,:) = f(w(1),w(2),w(3),t(1),t(2),t(3),X_wrd_inhomo(1,i),X_wrd_inhomo(2,i),X_wrd_inhomo(3,i));
end
end