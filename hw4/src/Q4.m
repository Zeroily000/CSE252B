clear;clc
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
%% DLT H
A = zeros(2*n,9);
for i = 1:n
    x1i = x1_norm(:,i);
    x2i = x2_norm(:,i);
    v = x2i + sign(x2i(1))*norm(x2i,2)*[1;0;0];
    Hv = eye(3) - 2 * (v*v')/(v'*v);
    m = Hv(2:3,:)';
    A(2*i-1:2*i,:) = kron(m',x1i');
end
[~,~,V] = svd(A);
h = V(:,end);
H_DLT = T2\reshape(h,3,3)'*T1;
H_DLT = H_DLT/norm(H_DLT,'fro');
save('../data/H_DLT.mat','H_DLT');

format longg
disp('H_DLT = ')
disp(H_DLT)
