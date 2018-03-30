function [A,B1,B2] = jcb(wu,wv,s,X_hat,f_A,f_B1,f_B2)
n = size(X_hat,2);
A = zeros(2,7,n);
B1 = zeros(2,3,n);
B2 = zeros(2,3,n);
for i = 1:n
    A(:,:,i) = f_A(X_hat(1,i),X_hat(2,i),X_hat(3,i),s,...
                   wu(1),wu(2),wu(3),wv(1),wv(2),wv(3));
    B1(:,:,i) = f_B1(X_hat(1,i),X_hat(2,i),X_hat(3,i));
    B2(:,:,i) = f_B2(X_hat(1,i),X_hat(2,i),X_hat(3,i),s,...
                     wu(1),wu(2),wu(3),wv(1),wv(2),wv(3));
end
end




