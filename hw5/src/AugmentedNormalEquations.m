function [delta_a,delta_b]=AugmentedNormalEquations(U,V,W,epsilon_a,epsilon_b,lambda)
n = size(V,3);
delta_b = zeros(size(epsilon_b,1),n);
S = U + lambda*eye(size(U));
e = epsilon_a;
for i = 1:n
    S = S-W(:,:,i)*((V(:,:,i)+lambda*eye(size(V(:,:,i))))\W(:,:,i)');
    e = e-W(:,:,i)*((V(:,:,i)+lambda*eye(size(V(:,:,i))))\epsilon_b(:,i));
end
delta_a = S\e;
for i = 1:n
    delta_b(:,i)=(V(:,:,i)+lambda*eye(size(V(:,:,i))))\...
                 (epsilon_b(:,i)-W(:,:,i)'*delta_a);
end
end