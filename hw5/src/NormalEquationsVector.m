function [epsilon_a,epsilon_b] = NormalEquationsVector(A,B1,B2,epsilon1,epsilon2)
n = size(A,3);
epsilon_a = zeros(size(A,2),size(epsilon2,2));
epsilon_b = zeros(size(B1,2),n);
for i = 1:n
    epsilon_a = epsilon_a + A(:,:,i)'*epsilon2(2*i-1:2*i);
    epsilon_b(:,i) = B1(:,:,i)'*epsilon1(2*i-1:2*i)+B2(:,:,i)'*epsilon2(2*i-1:2*i);
end
end