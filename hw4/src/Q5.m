clear;clc
load('../data/x_inlier.mat');
load('../data/H_DLT.mat');
n = size(x1_inlier_inhomo,2);
%% Jacobian
f_A = matlabFunction(jcbA);
f_B1 = matlabFunction(jcbB1);
f_B2 = matlabFunction(jcbB2);
%% Initialize parameter vector
% h_hat
h_bar = reshape(H_DLT',9,[]);
h_hat = parameterization(h_bar);
% x_hat
delta = SampsonError(x1_inlier_inhomo,x2_inlier_inhomo,H_DLT);
x_inhomo = x1_inlier_inhomo + delta(1:2,:);
x_bar = padarray(x_inhomo,[1 0],1,'post');
x_hat = parameterization(x_bar);
% sigma
sigma1 = eye(2*n);
sigma2 = eye(2*n);
% cost
previous_cost = inf;
% tolerance
tolerance = 0.0001;
%% LM
itr = 0;
disp('LM:')
fprintf('itr\tcost\n')
fprintf('-----------------\n')
% step_1
lambda = 0.001;
x1_hat = deparameterization(x_hat);
x2_hat = reshape(deparameterization(h_hat),3,3)' * deparameterization(x_hat);
x1_hat_inhomo = x1_hat(1:2,:) ./ x1_hat(3,:);
x2_hat_inhomo = x2_hat(1:2,:) ./ x2_hat(3,:);
epsilon1 = reshape(x1_inlier_inhomo - x1_hat_inhomo,[],1);
epsilon2 = reshape(x2_inlier_inhomo - x2_hat_inhomo,[],1);
current_cost = epsilon1'*(sigma1\epsilon1) + epsilon2'*(sigma2\epsilon2);
fprintf('%d\t%.4f\n', itr, current_cost)
% step_2
[A,B1,B2] = jcb(h_hat, x_hat, f_A, f_B1, f_B2);
while tolerance < previous_cost - current_cost || previous_cost < current_cost
    itr = itr+1;
    % step_3_4
    [U,V,W] = NormalEquationsMatrix(A,B1,B2);
    [epsilon_a,epsilon_b] = NormalEquationsVector(A,B1,B2,epsilon1,epsilon2);
    [delta_a,delta_b] = AugmentedNormalEquations(U,V,W,epsilon_a,epsilon_b,lambda);
    % step_5
    h_hat0 = h_hat + delta_a;
    x_hat0 = x_hat + delta_b;
    % step_6
    x1_hat0 = deparameterization(x_hat0);
    x2_hat0 = reshape(deparameterization(h_hat0),3,3)' * deparameterization(x_hat0);
    x1_hat_inhomo0 = x1_hat0(1:2,:) ./ x1_hat0(3,:);
    x2_hat_inhomo0 = x2_hat0(1:2,:) ./ x2_hat0(3,:);
    epsilon10 = reshape(x1_inlier_inhomo - x1_hat_inhomo0,[],1);
    epsilon20 = reshape(x2_inlier_inhomo - x2_hat_inhomo0,[],1);
    % step_7
    previous_cost = current_cost;
    current_cost = epsilon10'*(sigma1\epsilon10) + epsilon20'*(sigma2\epsilon20);
    fprintf('%d\t%.4f\n', itr, current_cost)
    if current_cost < previous_cost
        h_hat = h_hat0;
        x_hat = x_hat0;
        epsilon1 = epsilon10;
        epsilon2 = epsilon20;
        lambda = 0.1*lambda;
        [A,B1,B2] = jcb(h_hat, x_hat, f_A, f_B1, f_B2);
    else
        lambda = 10*lambda;
    end
end
fprintf('-----------------\n\n')
%% H_LM
H_LM = reshape(deparameterization(h_hat0),3,3)';
H_LM = H_LM/norm(H_LM,'fro');

format longg
disp('H_LM = ')
disp(H_LM)