clear;clc;
I0 = imread('../data/IMG_5030.JPG');
I1 = imread('../data/IMG_5031.JPG');
%% Corner Detection
fprintf('Detecting corners... ')
win_detect = 9;
win_spr = 9;
eigen_th = 5.57;
[r0,c0,x_f0,y_f0] = CornerCoordinate(I0,win_detect,win_spr,eigen_th);
[r1,c1,x_f1,y_f1] = CornerCoordinate(I1,win_detect,win_spr,eigen_th);
save('../data/c.mat','x_f0','y_f0','x_f1','y_f1');
fprintf('Done')
%% Plot
figure
imshow(I0);
% saveas(gcf,'Input0.png');
figure
imshow(I1);
% saveas(gcf,'Input1.png');
figure
imshow(I0);hold on
plot(y_f0,x_f0,'bs','MarkerSize',win_detect);hold off
% saveas(gcf,'CornerDetection0.png');
figure
imshow(I1);hold on
plot(y_f1,x_f1,'bs','MarkerSize',win_detect);hold off
% saveas(gcf,'CornerDetection1.png');
%% Print
disp(['Features in the first image: ',num2str(numel(x_f0))]);
disp(['Features in the second image: ',num2str(numel(x_f1))]);