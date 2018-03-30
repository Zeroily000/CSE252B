close all; clear; clc;
load('../data/x.mat');
I0 = imread('../data/IMG_5030.JPG');
I1 = imread('../data/IMG_5031.JPG');
x1 = padarray(x1_inhomo,[1 0],1,'post');
x2 = padarray(x2_inhomo,[1 0],1,'post');
n = size(x1,2);
%%
syms alpha 
syms a1 a2 a3 a4 a5 a6 a7 a8 a9
syms b1 b2 b3 b4 b5 b6 b7 b8 b9
tmp_F1 = [a1,a4,a7
          a2,a5,a8
          a3,a6,a9];
tmp_F2 = [b1,b4,b7
          b2,b5,b8
          b3,b6,b9];
t = det(alpha*tmp_F1+tmp_F2);
f = matlabFunction(flip(coeffs(t,alpha)));
%% mSAC
consensus_min_cost = inf;
max_trials = inf;
trials = 0;
threshold = 0;
tolerance = chi2inv(0.95,1) * 1;
p = 0.99;
s = 7;
rng(4)
while trials < max_trials && consensus_min_cost > threshold
    %% Select a random sample
    t = randperm(n,s)';
    x1_rand =  x1(:,t);
    x2_rand =  x2(:,t);
    % x1
    mu_x1 = mean(x1_rand,2);
    sigma_x1 = sum(var(x1_rand,0,2));
    s_x1 = sqrt(2/sigma_x1);
    T1 = [s_x1      0   -mu_x1(1)*s_x1
        0       s_x1   -mu_x1(2)*s_x1
        0        0        1      ];
    x1_norm = T1*x1_rand;
    % x2
    mu_x2 = mean(x2_rand,2);
    sigma_x2 = sum(var(x2_rand,0,2));
    s_x2 = sqrt(2/sigma_x2);
    T2 = [s_x2      0   -mu_x2(1)*s_x2
        0       s_x2   -mu_x2(2)*s_x2
        0        0        1      ];
    x2_norm = T2*x2_rand;
    %% Calculate model
    A = zeros(s,9);
    for i = 1:s
        A(i,:) = kron(x2_norm(:,i)', x1_norm(:,i)');
    end
    [~,~,V] = svd(A);
    F1 = reshape(V(:,end),  3,3)';
    F2 = reshape(V(:,end-1),3,3)';
    alpha = roots(f(F1(1),F1(2),F1(3),F1(4),F1(5),F1(6),F1(7),F1(8),F1(9),...
                    F2(1),F2(2),F2(3),F2(4),F2(5),F2(6),F2(7),F2(8),F2(9)));
    for i = 1:3
        if isreal(alpha(i))
            F_norm = alpha(i)*F1 + F2;
            break
        end
    end
    F = T2'*F_norm*T1;
    %% Error for each point
    error = SampsonError(x1,x2,F);
    %% Calculate cost
    cost = sum(error .* (error <= tolerance) + tolerance * (error > tolerance));
    %% Update maximum trials
    if cost < consensus_min_cost
        consensus_min_cost = cost;
        F_min_cost = F;
        inliers = error <= tolerance;
        w = sum(inliers) / n;
        max_trials = log(1-p) / log(1-w^s);
    end
    trials = trials + 1;
end
%% Save in- & out-liers
x1_inlier_inhomo = x1_inhomo(:,inliers);
x2_inlier_inhomo = x2_inhomo(:,inliers);
save('../data/x_inlier.mat','x1_inlier_inhomo','x2_inlier_inhomo');
x1_inhomo(:,inliers) = [];
x2_inhomo(:,inliers) = [];
x1_outlier_inhomo = x1_inhomo;
x2_outlier_inhomo = x2_inhomo;
save('../data/x_outlier.mat','x1_outlier_inhomo','x2_outlier_inhomo');
%% Plot
figure
imshow(insertShape(I0,'Line',[x1_inlier_inhomo',x2_inlier_inhomo'],...
    'Color','blue'));hold on
plot(x1_inlier_inhomo(1,:)',x1_inlier_inhomo(2,:)',...
    'bs','MarkerSize',win_match);hold off
% saveas(gcf,'FeatureMatch_inlier0.png');

figure
imshow(insertShape(I1,'Line',[x2_inlier_inhomo',x1_inlier_inhomo'],...
    'Color','blue'));hold on
plot(x2_inlier_inhomo(1,:)',x2_inlier_inhomo(2,:)',...
    'bs','MarkerSize',win_match);hold off
% saveas(gcf,'FeatureMatch_inlier1.png');
%% Print
clc;
disp(['The number of inliers is ',num2str(sum(inliers))]);
disp(['The value of max_trials is ',num2str(max_trials),]);
disp(['The the number of iteration is ',num2str(trials)]);