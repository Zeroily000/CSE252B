function [delta_a,delta_b] = AugmentedNormalEquations(U,V,W,epsilon_a,epsilon_b,lambda)
n = size(V,3);
delta_b = zeros(2,n);
S = U + lambda*eye(8);
e = epsilon_a;
for i = 1:n
    S = S - W(:,:,i)*((V(:,:,i) + lambda*eye(2))\W(:,:,i)');
    e = e - W(:,:,i)*((V(:,:,i) + lambda*eye(2))\epsilon_b(:,i));
end
delta_a = S\e;
for i = 1:n
    delta_b(:,i) = (V(:,:,i) + lambda*eye(2))\(epsilon_b(:,i) - W(:,:,i)'*delta_a);
end
end