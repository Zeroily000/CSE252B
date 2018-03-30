clear;clc
%% Load
load('../data/x_inlier.mat');
load('../data/F_DLT.mat');
n = size(x1_inlier_inhomo,2);
x1 = padarray(x1_inlier_inhomo,[1 0],1,'post');
x2 = padarray(x2_inlier_inhomo,[1 0],1,'post');
%% Initialize parameter vector
% wu,wv,s
[wu,wv,s] = Fparameterization(F_DLT);
% P
P1 = [eye(3),zeros(3,1)];
P2 = F2P(wu,wv,s);
% X_hat
g = polyg;
X = TwoViewTriangulation(x1,x2,wu,wv,s,g);
X_hat = Xparameterization(X);
%% LM
itr = 0;
disp('LM:')
fprintf('itr\tcost\n')
fprintf('-----------------\n')
% cost
previous_cost = inf;
% tolerance
tolerance = 0.0001;
% step_1
lambda = 0.001;
x1_hat = P1 * Xdeparameterization(X_hat);
x2_hat = P2 * Xdeparameterization(X_hat);
x1_hat_inhomo = x1_hat(1:2,:) ./ x1_hat(3,:);
x2_hat_inhomo = x2_hat(1:2,:) ./ x2_hat(3,:);
epsilon1 = reshape(x1_inlier_inhomo - x1_hat_inhomo,[],1);
epsilon2 = reshape(x2_inlier_inhomo - x2_hat_inhomo,[],1);
% step_2
[f_A,f_B1,f_B2] = jcbFunction;
[A,B1,B2] = jcb(wu,wv,s,X_hat,f_A,f_B1,f_B2);
% initial cost
init_cost = epsilon1'*epsilon1 + epsilon2'*epsilon2;
current_cost = init_cost;
fprintf('%d\t%.4f\n', itr, current_cost)
while tolerance < previous_cost - current_cost || previous_cost < current_cost
    % step_3_4
    [U,V,W] = NormalEquationsMatrix(A,B1,B2);
    [epsilon_a,epsilon_b] = NormalEquationsVector(A,B1,B2,epsilon1,epsilon2);
    [delta_a,delta_b] = AugmentedNormalEquations(U,V,W,epsilon_a,epsilon_b,lambda);
    % step_5
    wu0 = wu + delta_a(1:3);
    wv0 = wv + delta_a(4:6);
    s0 = s + delta_a(7);
    X_hat0 =  X_hat + delta_b;
    % step_6
    x1_hat0 = P1 * Xdeparameterization(X_hat0);
    x2_hat0 = F2P(wu0,wv0,s0) * Xdeparameterization(X_hat0);       
    x1_hat_inhomo0 = x1_hat0(1:2,:) ./ x1_hat0(3,:);
    x2_hat_inhomo0 = x2_hat0(1:2,:) ./ x2_hat0(3,:);
    epsilon10 = reshape(x1_inlier_inhomo - x1_hat_inhomo0,[],1);
    epsilon20 = reshape(x2_inlier_inhomo - x2_hat_inhomo0,[],1);
    % step_7
    previous_cost = current_cost;
    current_cost = epsilon10'*epsilon10 + epsilon20'*epsilon20;
    if current_cost < previous_cost
        wu = wu0;
        wv = wv0;
        s = s0;
        X_hat = X_hat0;
        epsilon1 = epsilon10;
        epsilon2 = epsilon20;
        lambda = 0.1*lambda;
        [A,B1,B2] = jcb(wu,wv,s,X_hat,f_A,f_B1,f_B2);
        if current_cost < init_cost
            itr = itr + 1;
            fprintf('%d\t%.4f\n', itr, current_cost)
        end
    else
        lambda = 10*lambda;
    end
end
fprintf('-----------------\n\n')
%% F_LM
F_LM = Fdeparameterization(wu,wv,s);
F_LM = F_LM/norm(F_LM,'fro');
save('../data/F_LM.mat','F_LM');

format longg
disp('F_LM = ')
disp(F_LM)
