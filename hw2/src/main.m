close all;
clear;
clc;
%% Load Data
x = load('hw2_points2D.txt');
X = load('hw2_points3D.txt');
n = numel(x);
A = zeros(n,12);
%% Data Normalization
mu_x = mean(x);
sigma_x = sum(var(x));
s_x = sqrt(2/sigma_x);
T = [s_x      0   -mu_x(1)*s_x
     0       s_x  -mu_x(2)*s_x
     0        0        1      ];
 
mu_X = mean(X);
sigma_X = sum(var(X));
s_X = sqrt(3/sigma_X);
U = [s_X      0        0     -mu_X(1)*s_X
     0       s_X       0     -mu_X(2)*s_X
     0        0       s_X    -mu_X(3)*s_X
     0        0        0         1     ];
x_bar = padarray(x,[0 1],1,'post');
X_bar = padarray(X,[0 1],1,'post');
x_bar_norm = T*x_bar';
X_bar_norm = U*X_bar';
%% P
for i = 1:2:n
    xi = x_bar_norm(:,(i+1)/2);
    Xi = X_bar_norm(:,(i+1)/2);
    v = xi + sign(xi(1))*norm(xi,2)*[1;0;0];
    Hv = eye(3) - 2 * (v*v')/(v'*v);
    m = Hv(2:3,:)';
    A(i:i+1,:) = kron(m',Xi');
end
[~,~,V] = svd(A);
p = V(:,end);
format shortg
P_DLT = T\reshape(p,4,3)'*U;
P_DLT = P_DLT/norm(P_DLT,'fro');
disp('P_DLT = ')
disp(P_DLT)
%% Vecterization
x_norm = reshape(x_bar_norm(1:2,:),[],1);
%% Initialization
% P_init = reshape(p,4,3)';
% p_hat
p_bar = p;
p_hat = parameterization(p_bar);
% x_hat
x_hat = reshape(p,4,3)'*X_bar_norm;
x_hat = x_hat(1:2,:) ./ x_hat(3,:);
x_hat = reshape(x_hat,[],1);
% sigma
sigma = T(1,1)^2 * eye(n);
% cost
previous_cost = inf;
tolerance = 0.0001;
%% LM
i = 0;
disp('LM:')
fprintf('itr\tcost\n')
fprintf('-----------------\n')
% step_1
lambda = 0.001;
epsilon = x_norm - x_hat;
format short
current_cost = epsilon'*(sigma\epsilon);
fprintf('%d\t%.4f\n', i, current_cost)
% step_2
J = jcb(x_hat,X_bar_norm,p_hat);
while tolerance < previous_cost - current_cost
    i = i + 1;
    % step_3_4
    delta = (J'*(sigma\J)+lambda*eye(11))\(J'*(sigma\epsilon));
    % step_5
    p_hat0 = p_hat + delta;
    % step_6
    x_hat0 = estimate(p_hat0,X_bar_norm);
    epsilon0 = x_norm - x_hat0;
    % step_7
    format short
    previous_cost = epsilon'*(sigma\epsilon);
    current_cost = epsilon0'*(sigma\epsilon0);
    fprintf('%d\t%.4f\n', i, current_cost)
    if current_cost < previous_cost
        p_hat = p_hat0;
        epsilon = epsilon0;
        lambda = 0.1*lambda;
        J = jcb(x_hat0,X_bar_norm,p_hat);
    else
        lambda = 10*lambda;
    end
end

p_bar_final = deparameterization(p_hat0);
format shortg
P_LM = T\reshape(p_bar_final,4,3)'*U;
P_LM = P_LM/norm(P_LM,'fro');
fprintf('\nP_LM = \n')
disp(P_LM)