function [r,c,x_f1,x_f2] = CornerCoordinate(im,win1,win2,threshold)
im = double(rgb2gray(im));
A = zeros(size(im));
k = [-1 8 0 -8 1]' / 12;
imx = imfilter(im,k','symmetric');
imy = imfilter(im,k,'symmetric');
%% Gradient matrix
for i = (win1+1)/2 : size(im,1) - (win1-1)/2
    for j = (win1+1)/2 : size(im,2) - (win1-1)/2
        im_win_x = imx(i - (win1-1)/2:i + (win1-1)/2,...
                       j - (win1-1)/2:j + (win1-1)/2);
        im_win_y = imy(i - (win1-1)/2:i + (win1-1)/2,...
                       j + (win1-1)/2:j + (win1-1)/2);

        N = zeros(2,2);
        N(1,1) = sum(sum(im_win_x.^2));
        N(1,2) = sum(sum(im_win_x .* im_win_y));
        N(2,1) = N(1,2);
        N(2,2) = sum(sum(im_win_y.^2));
        N = N/win1^2;
        lambda = (trace(N) - sqrt(trace(N)^2 - 4*det(N)))/2;
        if lambda > threshold
            A(i,j) = lambda;
        end
    end
end
%% Non-maximum suppression 
B = ordfilt2(A,win1^2,ones(win1,win1));
C = (B == A & B ~= 0);
[r,c] = find(C);
x_f1 = zeros(size(r));
x_f2 = zeros(size(r));
%% find corner
for k = 1:numel(r)
    T = zeros(2,2);
    y = zeros(2,1);
    im_win_x = imx(r(k) - (win2-1)/2:r(k) + (win2-1)/2,...
                   c(k) - (win2-1)/2:c(k) + (win2-1)/2);
    im_win_y = imy(r(k) - (win2-1)/2:r(k) + (win2-1)/2,...
                   c(k) - (win2-1)/2:c(k) + (win2-1)/2);
    xx = (r(k) - (win2-1)/2 : r(k) + (win2-1)/2)' .* ones(win2);
    yy = (c(k) - (win2-1)/2 : c(k) + (win2-1)/2) .* ones(win2);
    
    y(1) = sum(sum( xx.*im_win_x.^2 + yy.*im_win_x.*im_win_y ));
    y(2) = sum(sum( xx.*im_win_x.*im_win_y + yy.*im_win_y.^2 ));
    T(1,1) = sum(sum(im_win_x.^2));
    T(1,2) = sum(sum(im_win_x .* im_win_y));
    T(2,1) = T(1,2);
    T(2,2) = sum(sum(im_win_y.^2));
    x = T\y;
    
    x_f1(k) = x(1);
    x_f2(k) = x(2);
end
end