function J = jcb(x_hat,X_bar_norm,p_hat)
% inputs:
% 1. inhomogeneous 2D points x vector;
% 2. homogeneous 3D points X;
% 3. p_bar(12*1)
p_bar = deparameterization(p_hat);
p_hat_norm = norm(p_hat,2);
r = size(x_hat,1);
c = numel(p_hat);
J = zeros(r,c);
for i = 1:2:r
    w = p_bar(end-3:end)'*X_bar_norm(:,(i+1)/2);
    dx_dpbar = [X_bar_norm(:,(i+1)/2)' zeros(1,4) -x_hat(i)*X_bar_norm(:,(i+1)/2)' % 2*12
                zeros(1,4) X_bar_norm(:,(i+1)/2)' -x_hat(i+1)*X_bar_norm(:,(i+1)/2)']/w;
    da_dphat = -p_bar(2:end)'/2;
    if p_hat_norm == 0
        db_dphat = eye(c)/2;
    else
        db_dphat = sin(p_hat_norm/2)/p_hat_norm*eye(c) + ...
            (p_hat_norm/2*cos(p_hat_norm/2)-sin(p_hat_norm/2))/p_hat_norm^3*(p_hat*p_hat');
    end
    dpbar_dphat = [da_dphat;db_dphat];
    J(i:i+1,:) = dx_dpbar*dpbar_dphat;
end