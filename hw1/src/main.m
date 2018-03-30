close all;
clear;
clc;
%% Corner Detection
fprintf('Detecting corners... ')
win_detect = 11;
win_spr = 9;
eigen_th = 5;
I0 = imread('../data/price_center20.JPG');
I1 = imread('../data/price_center21.JPG');
[r0,c0,x_f0,y_f0] = CornerCoordinate(I0,win_detect,win_spr,eigen_th);
[r1,c1,x_f1,y_f1] = CornerCoordinate(I1,win_detect,win_spr,eigen_th);
fprintf('Done\n')

%% Feature Matching
fprintf('Matching features... ')
win_match = 11;
simi_th = 0.6;
dist_th = 0.9;
[X0,Y0,X1,Y1] = FeatureMatch(I0,I1,x_f0,y_f0,x_f1,y_f1,win_match,simi_th,dist_th);
fprintf('Done\n')

%% Plots
figure(1)
subplot(1,2,1)
imshow(I0);
xlabel('(a)');
title('''price\_center20''');
subplot(1,2,2)
imshow(I1);
xlabel('(b)');
title('''price\_center21''');

figure(2)
subplot(1,2,1)
imshow(I0);hold on
plot(c0,r0,'bs','MarkerSize',win_detect);hold off
xlabel('(a)');
title('Corner Detection of ''price\_center20''');
subplot(1,2,2)
imshow(I1);hold on
xlabel('(b)')
plot(c1,r1,'bs','MarkerSize',win_detect);hold off
title('Corner Detection of ''price\_center21''');

figure(3)
subplot(1,2,1)
imshow(I0);hold on
plot(y_f0,x_f0,'bs','MarkerSize',win_detect);hold off
xlabel('(a)');
title('Corner Detection of ''price\_center20''');
subplot(1,2,2)
imshow(I1);hold on
xlabel('(b)')
plot(y_f1,x_f1,'bs','MarkerSize',win_detect);hold off
title('Corner Detection of ''price\_center21''');

figure(4)
subplot(1,2,1)
imshow(insertShape(I0,'Line',[Y0 X0 Y1 X1],'Color','blue'));
hold on
plot(Y0,X0,'bs','MarkerSize',win_match);
hold off
xlabel('(a)');
title('Feature Matching of ''price\_center20''');
subplot(1,2,2)
imshow(insertShape(I1,'Line',[Y1 X1 Y0 X0],'Color','blue'));
hold on
plot(Y1,X1,'bs','MarkerSize',win_match);
hold off
xlabel('(b)');
title('Feature Matching of ''price\_center21''');