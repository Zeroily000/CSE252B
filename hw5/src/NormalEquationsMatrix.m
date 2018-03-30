function [U,V,W] = NormalEquationsMatrix(A,B1,B2)
n = size(A,3);
U = zeros(size(A,2),size(A,2));
V = zeros(size(B1,2),size(B1,2),n);
W = zeros(size(A,2),size(B1,2),n);
for i = 1:n
    U = U + A(:,:,i)' * A(:,:,i);
    V(:,:,i) = B1(:,:,i)'*B1(:,:,i) + B2(:,:,i)'*B2(:,:,i);
    W(:,:,i) = A(:,:,i)'*B2(:,:,i);
end
end