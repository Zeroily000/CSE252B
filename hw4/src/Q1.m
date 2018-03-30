close all; clear; clc;
I0 = imread('../data/price_center20.JPG');
I1 = imread('../data/price_center21.JPG');
%% Corner Detection
fprintf('Detecting corners... ')
win_detect = 9;
win_spr = 9;
eigen_th = 6.5;
[r0,c0,x_f0,y_f0] = CornerCoordinate(I0,win_detect,win_spr,eigen_th);
[r1,c1,x_f1,y_f1] = CornerCoordinate(I1,win_detect,win_spr,eigen_th);
save('../data/c.mat','x_f0','y_f0','x_f1','y_f1');
fprintf('Done\n')
disp(['The number of feature detected in ''price\_center20'' is '...
    ,num2str(numel(r0))]);
disp(['The number of feature detected in ''price\_center21'' is '...
    ,num2str(numel(r1))]);
%% Plot
figure
imshow(I0);
title('price\_center20');
% saveas(gcf,'PriceCenter_20.png');
figure
imshow(I1);
title('price\_center21');
% saveas(gcf,'PriceCenter_21.png');
figure
imshow(I0);hold on
plot(y_f0,x_f0,'bs','MarkerSize',win_detect);hold off
title('Corner Detection of ''price\_center20''');
% saveas(gcf,'CornerDetection_20.png');
figure
imshow(I1);hold on
plot(y_f1,x_f1,'bs','MarkerSize',win_detect);hold off
title('Corner Detection of ''price\_center21''');
% saveas(gcf,'CornerDetection_21.png');