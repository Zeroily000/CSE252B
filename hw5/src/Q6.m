close all; clear; clc;
%% Load Data
load('../data/x_outlier.mat');
load('../data/F_LM.mat');
I0 = imread('../data/IMG_5030.JPG');
I1 = imread('../data/IMG_5031.JPG');
x1 = padarray(x1_outlier_inhomo,[1 0],1,'post');
n = size(x1,2);
%% Mapping
rng(1)
t = randperm(n,3)';
x1_rnd =  x1(:,t);
x2_rnd_inhomo = x2_outlier_inhomo(:,t);
l2 = F_LM*x1_rnd;
%% Plot
figure
imshow(rgb2gray(I0));hold on
plot(x1_outlier_inhomo(1,t)',x1_outlier_inhomo(2,t)','ys','MarkerSize',7);
hold off
% saveas(gcf,'Outlier0.png');
figure
x = (1:size(I0,2))';
y = (-l2(1,:) .* x - l2(3,:)) ./ l2(2,:);
imshow(rgb2gray(I1));hold on
plot(x2_rnd_inhomo(1,:)',x2_rnd_inhomo(2,:)','ys','MarkerSize',7);
plot(x,y,'y','LineWidth',2);hold off
% saveas(gcf,'Outlier1.png');
%% Distance
clc;
distance=abs(l2(1,:).*x2_rnd_inhomo(1,:)+l2(2,:).*x2_rnd_inhomo(2,:)+l2(3,:))...
         ./ sqrt(l2(1,:).^2+l2(2,:).^2);
disp(['The distance between the first point and the line is ',...
    num2str(distance(1))]);
disp(['The distance between the second point and the line is ',...
    num2str(distance(2))]);
disp(['The distance between the third point and the line is ',...
    num2str(distance(3))]);