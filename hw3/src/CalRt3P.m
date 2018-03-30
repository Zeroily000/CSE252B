function P_hat = CalRt3P(X1_wld,X2_wld,X3_wld,X1_cam,X2_cam,X3_cam)
P_hat = [];
n = 3;
[m,soln] = size(X1_cam);

X = [X1_wld,X2_wld,X3_wld];
mu_x = mean(X,2);
sigma_x = norm(X-mu_x,'fro')^2/n;
err = inf;
for i = 1:soln
    Y = [X1_cam(:,i),X2_cam(:,i),X3_cam(:,i)];
    mu_y = mean(Y,2);
    Sigma_xy = (Y-mu_y)*(X-mu_x)'/n;
    [U,D,V] = svd(Sigma_xy);
    S = eye(m);  
    if rank(Sigma_xy) < m-1
        continue
    elseif rank(Sigma_xy) == m-1
        if round(det(U)*det(V)) == -1
            S(end) = -1;
        end
    else
        if det(Sigma_xy) < 0
            S(end) = -1;
        end
    end
    R = U*S*V';
    c = trace(D*S)/sigma_x;
    t = mu_y - c*R*mu_x;
    if norm(Y - c*R*X - t)^2 < err
        err = norm(Y - c*R*X - t)^2;
        P_hat = [R,t];
    end
end
end