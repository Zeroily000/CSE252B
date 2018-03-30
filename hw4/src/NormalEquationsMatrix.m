function [U,V,W] = NormalEquationsMatrix(A,B1,B2)
n = size(A,3);
U = zeros(8,8);
V = zeros(2,2,n);
W = zeros(8,2,n);
for i = 1:n
    U = U + A(:,:,i)' * A(:,:,i);
    V(:,:,i) = B1(:,:,i)'*B1(:,:,i) + B2(:,:,i)'*B2(:,:,i);
    W(:,:,i) = A(:,:,i)'*B2(:,:,i);
end
end