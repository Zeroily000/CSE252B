close all; clear; clc;
load('../data/c.mat');
I0 = imread('../data/price_center20.JPG');
I1 = imread('../data/price_center21.JPG');
%% Feature Matching
fprintf('Matching features... ')
win_match = 11;
simi_th = 0.7;
dist_th = 0.9;
[X0,Y0,X1,Y1] = FeatureMatch(I0,I1,x_f0,y_f0,x_f1,y_f1,...
                             win_match,simi_th,dist_th);
fprintf('Done\n')
disp(['The number of matched features is ',num2str(numel(X0))]);
x1_inhomo = [Y0';X0'];
x2_inhomo = [Y1';X1'];
save('../data/x.mat','x1_inhomo','x2_inhomo','win_match');
%% Plot
figure
imshow(insertShape(I0,'Line',[Y0 X0 Y1 X1],'Color','blue'));hold on
plot(Y0,X0,'bs','MarkerSize',win_match);hold off
title('Feature Matching of ''price\_center20''');
% saveas(gcf,'FeatureMatch_20.png');
figure
imshow(insertShape(I1,'Line',[Y1 X1 Y0 X0],'Color','blue'));hold on
plot(Y1,X1,'bs','MarkerSize',win_match);hold off
title('Feature Matching of ''price\_center21''');
% saveas(gcf,'FeatureMatch_21.png');