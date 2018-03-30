close all; clear; clc;
%% Load Data
x_img_inhomo = load('../data/hw3_points2D.txt')';
X_wld_inhomo = load('../data/hw3_points3D.txt')';
x_img_homo = padarray(x_img_inhomo,[1,0],1,'post');
X_wld_homo = padarray(X_wld_inhomo,[1,0],1,'post');
K = [1545.0966799187809,0,639.5;0,1545.0966799187809,359.5;0,0,1];
x_img_norm_homo = K\x_img_homo;
x_img_norm_homo = x_img_norm_homo./sign(x_img_norm_homo(3,:)) / norm(x_img_norm_homo);
x_img_norm_inhomo = x_img_norm_homo(1:2,:) ./ x_img_norm_homo(3,:);
n = size(x_img_inhomo,2);
%% mSAC
rng(53)
consensus_min_cost = inf;
max_trials = inf;
trials = 0;
threshold = 0;
tolerance = chi2inv(0.95,2) * 1;
p = 0.99;
s = 3;
k = 1;
while trials < max_trials && consensus_min_cost > threshold
    %% Select a random sample
    i = randperm(n,3)';
    p1_wld = X_wld_inhomo(:,i(1));
    p2_wld = X_wld_inhomo(:,i(2));
    p3_wld = X_wld_inhomo(:,i(3));
    
    q1 = [x_img_norm_inhomo(:,i(1));1] ;
    q2 = [x_img_norm_inhomo(:,i(2));1] ;
    q3 = [x_img_norm_inhomo(:,i(3));1] ;
    
    j1 = q1/norm(q1);
    j2 = q2/norm(q2);
    j3 = q3/norm(q3);
    
    a = norm(p2_wld-p3_wld);
    b = norm(p1_wld-p3_wld);
    c = norm(p1_wld-p2_wld);
    %% Calculate model
    %         Finsterwalder
    [X1_cam,X2_cam,X3_cam]=Finsterwalder(a,b,c,j1,j2,j3);
    %         Projection Matrix
    P_hat = CalRt3P(p1_wld,p2_wld,p3_wld,X1_cam,X2_cam,X3_cam);
    %         Check if solution exist
    if ~numel(P_hat)
        continue
    end
    trials = trials + 1;
    %% Error for each point
    x_pro_homo = K * P_hat * X_wld_homo;
    x_pro_inhomo = x_pro_homo(1:2,:) ./ x_pro_homo(3,:);
    error = sum((x_img_inhomo - x_pro_inhomo).^2);
    %% Calculate cost
    cost = sum(error .* (error <= tolerance) + tolerance * (error > tolerance));
    %% Update maximum trials
    if cost < consensus_min_cost
        consensus_min_cost = cost;
        P_min_cost = P_hat;
        inliers = error <= tolerance;
        w = sum(inliers) / n;
        max_trials = log(1-p) / log(1-w^s);
    end
    k = k + 1;
end
fprintf('inliers:\t %d\n', sum(inliers))
fprintf('maximum trials:\t %.4f\n\n', max_trials)
%% b
x_img_norm_inlier_inhomo = x_img_norm_inhomo(:,inliers);
X_wld_inlier_inhomo = X_wld_inhomo(:,inliers);
X_wld_inlier_homo = X_wld_homo(:,inliers);
n = size(x_img_norm_inlier_inhomo,2);
%% Control points in world coordinate frame
mu_X_wld = mean(X_wld_inlier_inhomo,2);
Sigma_X_wld = cov(X_wld_inlier_inhomo');
[~,~,V] = svd(Sigma_X_wld);
var_X_wld = trace(Sigma_X_wld);
s = sqrt(var_X_wld / 3);
C1_wld_inhomo = mu_X_wld;
C2_wld_inhomo = s*V(:,1) + mu_X_wld;
C3_wld_inhomo = s*V(:,2) + mu_X_wld;
C4_wld_inhomo = s*V(:,3) + mu_X_wld;
%% Parameterize 3D points
A = [C2_wld_inhomo - C1_wld_inhomo,...
     C3_wld_inhomo - C1_wld_inhomo,...
     C4_wld_inhomo - C1_wld_inhomo];
b = X_wld_inlier_inhomo - C1_wld_inhomo;
X_prm_wrd_homo = [1-sum(A\b);A\b];
%% Control points in camera coordinate frame
m = zeros(2*n,12);
for i = 1:n
    m(2*i - 1 : 2*i,:) = [X_prm_wrd_homo(1,i) 0 -X_prm_wrd_homo(1,i)*x_img_norm_inlier_inhomo(1,i), X_prm_wrd_homo(2,i) 0 -X_prm_wrd_homo(2,i)*x_img_norm_inlier_inhomo(1,i), X_prm_wrd_homo(3,i) 0 -X_prm_wrd_homo(3,i)*x_img_norm_inlier_inhomo(1,i), X_prm_wrd_homo(4,i) 0 -X_prm_wrd_homo(4,i)*x_img_norm_inlier_inhomo(1,i)
                          0 X_prm_wrd_homo(1,i) -X_prm_wrd_homo(1,i)*x_img_norm_inlier_inhomo(2,i), 0 X_prm_wrd_homo(2,i) -X_prm_wrd_homo(2,i)*x_img_norm_inlier_inhomo(2,i), 0 X_prm_wrd_homo(3,i) -X_prm_wrd_homo(3,i)*x_img_norm_inlier_inhomo(2,i), 0 X_prm_wrd_homo(4,i) -X_prm_wrd_homo(4,i)*x_img_norm_inlier_inhomo(2,i)];
end
[~,~,V] = svd(m);
C1_cam_inhomo = V(1:3,end);
C2_cam_inhomo = V(4:6,end);
C3_cam_inhomo = V(7:9,end);
C4_cam_inhomo = V(10:12,end);
%% Deparameterize 3D points in camera coordinate frame
X_cam_inhomo = [C1_cam_inhomo,C2_cam_inhomo,C3_cam_inhomo,C4_cam_inhomo] * X_prm_wrd_homo;
mu_X_cam = mean(X_cam_inhomo,2);
Sigma_X_cam = cov(X_cam_inhomo');
var_X_cam = trace(Sigma_X_cam);
if mu_X_cam(3) < 0
    beta = -sqrt(var_X_wld/var_X_cam);
else
    beta = sqrt(var_X_wld/var_X_cam);
end
X_cam_inhomo = beta*X_cam_inhomo;
P_lin = CalRtnP(X_wld_inlier_inhomo,X_cam_inhomo);
R_lin = P_lin(:,1:3);
t_lin = P_lin(:,end);
format longg
disp('R_lin = ')
disp(R_lin)
disp('t_lin = ')
disp(t_lin)

%% LM
itr = 0;
disp('LM:')
fprintf('itr\tcost\n')
fprintf('-----------------\n')
%% Jacobian
syms w1 w2 w3 t1 t2 t3 X1 X2 X3
ww = [w1;w2;w3];
XX = [X1;X2;X3];
tt = [t1;t2;t3];
theta = norm(ww);
Xrotate_large = XX + sinc(theta/pi)*cross(ww,XX) + (1-cos(theta))/theta^2*cross(ww,cross(ww,XX));
x_homo_large = Xrotate_large + tt;
x_inhomo_large = x_homo_large(1:2)/x_homo_large(3);
f_large([w1 w2 w3 t1 t2 t3 X1 X2 X3]) = jacobian(x_inhomo_large,[ww.',tt.']);

Xrotate_small = XX + cross(ww,XX);
x_homo_small = Xrotate_small + tt;
x_inhomo_small = x_homo_small(1:2)/x_homo_small(3);
f_small([w1 w2 w3 t1 t2 t3 X1 X2 X3]) = jacobian(x_inhomo_small,[ww.',tt.']);
%% Initialization
x_prj_norm_inlier_homo = [R_lin,t_lin]*X_wld_inlier_homo;
x_prj_norm_inlier_inhomo = x_prj_norm_inlier_homo(1:2,:) ./ x_prj_norm_inlier_homo(3,:);

w_lin = parameterization(R_lin);
p_hat = [w_lin;t_lin];

previous_cost = inf;
tolerance = 0.00001;

% step_1
lambda = 0.001;
K_inv = inv(K);
sigma = eye(2*n);
for i = 1:n 
    sigma(2*i-1:2*i,2*i-1:2*i) = K_inv(1:2,1:2) * eye(2) * K_inv(1:2,1:2)';
end
epsilon = reshape(x_img_norm_inlier_inhomo - x_prj_norm_inlier_inhomo,[],1);




% format short
% iteration = 0;
init_cost = epsilon'*(sigma\epsilon);
current_cost = init_cost;
fprintf('%d\t%.4f\n', itr, current_cost)

J = jcb(f_large,f_small,X_wld_inlier_inhomo,w_lin,t_lin);
while tolerance < previous_cost - current_cost
    itr = itr + 1;
    % step_3_4
    delta = (J'*(sigma\J)+lambda*eye(6))\(J'*(sigma\epsilon));
    % step_5
    p_hat0 = p_hat + delta;
    w0 = p_hat0(1:3);
    t0 = p_hat0(end-2:end);
    % step_6
    R0 = deparameterization(w0);
    x_prj_norm_inlier_homo = [R0,t0]*X_wld_inlier_homo;
    x0_prj_norm_inlier_inhomo = x_prj_norm_inlier_homo(1:2,:) ./ x_prj_norm_inlier_homo(3,:);
    epsilon0 = reshape(x_img_norm_inlier_inhomo - x0_prj_norm_inlier_inhomo,[],1);
    % step_7
    format short
    previous_cost = epsilon'*(sigma\epsilon);
    current_cost = epsilon0'*(sigma\epsilon0);
    fprintf('%d\t%.4f\n', itr, current_cost)
    if current_cost < previous_cost
        p_hat = p_hat0;
        epsilon = epsilon0;
        lambda = 0.1*lambda;
        J = jcb(f_large,f_small,X_wld_inlier_inhomo,w0,t0);
    else
        lambda = 10*lambda;
    end
end
fprintf('-----------------\n\n')
w_LM = p_hat(1:3);
R_LM = deparameterization(w_LM);
t_LM = p_hat(end-2:end);
format longg
disp('w_LM = ')
disp(w_LM)
disp('R_LM = ')
disp(R_LM)
disp('t_LM = ')
disp(t_LM)
