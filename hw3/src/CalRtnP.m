function P_hat = CalRtnP(X_wld_inhomo,X_cam_inhomo)
[m,n] = size(X_cam_inhomo);

mu_x = mean(X_wld_inhomo,2);
mu_y = mean(X_cam_inhomo,2);
Sigma_xy = (X_cam_inhomo-mu_y)*(X_wld_inhomo-mu_x)'/n;
[U,~,V] = svd(Sigma_xy);
S = eye(m);
if rank(Sigma_xy) >= m-1
    if det(Sigma_xy) < 0
        S(end) = -1;
    end
    R = U*S*V';
    c = 1;
    t = mu_y - c*R*mu_x;
    P_hat = [R,t];
else
    P_hat = [];
end
end