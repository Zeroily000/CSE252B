clear;clc
%% Load
load('../data/x_inlier.mat');
x1 = padarray(x1_inlier_inhomo,[1 0],1,'post');
x2 = padarray(x2_inlier_inhomo,[1 0],1,'post');
n = size(x1,2);
%% Data Normalization
% x1
mu_x1 = mean(x1,2);
sigma_x1 = sum(var(x1,0,2));
s_x1 = sqrt(2/sigma_x1);
T1 = [s_x1      0   -mu_x1(1)*s_x1
    0       s_x1   -mu_x1(2)*s_x1
    0        0        1      ];
x1_norm = T1*x1;
% x2
mu_x2 = mean(x2,2);
sigma_x2 = sum(var(x2,0,2));
s_x2 = sqrt(2/sigma_x2);
T2 = [s_x2      0   -mu_x2(1)*s_x2
    0       s_x2   -mu_x2(2)*s_x2
    0        0        1      ];
x2_norm = T2*x2;
%% DLT F
A = zeros(n,9);
for i = 1:n
    A(i,:) = kron(x2_norm(:,i)', x1_norm(:,i)');
end
[~,~,V] = svd(A);
F_norm = reshape(V(:,end),3,3)';
%% Rank 2 constraint
[U,D,V] = svd(F_norm);
F_norm = U*diag([D(1,1),D(2,2),0])*V';
F_DLT = T2'*F_norm*T1;
%% Print
F_DLT = F_DLT/norm(F_DLT,'fro');
save('../data/F_DLT.mat','F_DLT');
format longg
disp('F_DLT = ')
disp(F_DLT)