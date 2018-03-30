clear;clc;
load('../data/x.mat');
x1 = padarray(x1_inhomo,[1 0],1,'post');
x2 = padarray(x2_inhomo,[1 0],1,'post');
n = size(x1,2);
%% mSAC
consensus_min_cost = inf;
max_trials = inf;
trials = 0;
threshold = 0;
tolerance = chi2inv(0.95,2) * 1;
p = 0.99;
s = 3;
k = 1;
rng(2)
while trials < max_trials && consensus_min_cost > threshold
    %% Select a random sample
    i = randperm(n,4)';
    x11 = x1(:,i(1));
    x12 = x1(:,i(2));
    x13 = x1(:,i(3));
    x14 = x1(:,i(4));
    
    x21 = x2(:,i(1));
    x22 = x2(:,i(2));
    x23 = x2(:,i(3));
    x24 = x2(:,i(4));
    %% Calculate model
    lambda1 = [x11,x12,x13] \ x14;
    H1 = inv([lambda1(1)*x11,lambda1(2)*x12,lambda1(3)*x13]);    
    lambda2 = [x21,x22,x23] \ x24;
    H2 = inv([lambda2(1)*x21,lambda2(2)*x22,lambda2(3)*x23]);   
    H = H2\H1;
    %% Error for each point
    delta = SampsonError(x1_inhomo,x2_inhomo,H);
    error = sum(delta.^2);
    %% Calculate cost
    cost = sum(error .* (error<=tolerance) + tolerance * (error>tolerance));
    %% Update maximum trials
    if cost < consensus_min_cost
        consensus_min_cost = cost;
        H_min_cost = H;
        inliers = error <= tolerance;
        w = sum(inliers) / n;
        max_trials = log(1-p) / log(1-w^s);
    end
    trials = trials + 1;
end
disp(['The number of inliers is ',num2str(sum(inliers))]);
disp(['The value of max_trials is ',num2str(max_trials),...
    ', so the the number of attempts is ',num2str(ceil(max_trials))]);
x1_inlier_inhomo = x1_inhomo(:,inliers);
x2_inlier_inhomo = x2_inhomo(:,inliers);
save('../data/x_inlier.mat','x1_inlier_inhomo','x2_inlier_inhomo');
%% Plot
I0 = imread('../data/price_center20.JPG');
I1 = imread('../data/price_center21.JPG');
figure
imshow(insertShape(I0,'Line',[x1_inlier_inhomo',x2_inlier_inhomo'],...
    'Color','blue'));hold on
plot(x1_inlier_inhomo(1,:)',x1_inlier_inhomo(2,:)',...
    'bs','MarkerSize',win_match);hold off
title('Feature Matching of ''price\_center20''');
% saveas(gcf,'FeatureMatch_20_inlier.png');

figure
imshow(insertShape(I1,'Line',[x2_inlier_inhomo',x1_inlier_inhomo'],...
    'Color','blue'));hold on
plot(x2_inlier_inhomo(1,:)',x2_inlier_inhomo(2,:)',...
    'bs','MarkerSize',win_match);hold off
title('Feature Matching of ''price\_center21''');
% saveas(gcf,'FeatureMatch_21_inlier.png');