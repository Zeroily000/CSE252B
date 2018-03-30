function [X0,Y0,X1,Y1] = FeatureMatch(im0,im1,x0,y0,x1,y1,win,simi_th,dist_th)
im0 = double(rgb2gray(im0));
im1 = double(rgb2gray(im1));
half_win = (win-1)/2;
%%
% xx0 = round(x0)+half_win;
% xx1 = round(x1)+half_win;
% yy0 = round(y0)+half_win;
% yy1 = round(y1)+half_win;
xx0 = x0+half_win;
xx1 = x1+half_win;
yy0 = y0+half_win;
yy1 = y1+half_win;
im0 = padarray(im0,[half_win,half_win],'symmetric');
im1 = padarray(im1,[half_win,half_win],'symmetric');
%% Xcorrelation matrix
xc = zeros(numel(xx0),numel(xx1));
crd = zeros(numel(xx0),numel(xx1));
for i = 1:numel(xx0)
    i
    X = fix(xx0(i)) - half_win:ceil(xx0(i)) + half_win;
    Y = fix(yy0(i)) - half_win:ceil(yy0(i)) + half_win;
    im0_win = im0(X,Y);
    [X,Y] = meshgrid(X,Y);
    [Xq,Yq] = meshgrid(xx0(i)-half_win:xx0(i)+half_win,...
                       yy0(i) - half_win:yy0(i) + half_win);
    im0_win_intep = interp2(X,Y,im0_win,Xq,Yq,'linear');    
%     im0_win = im0(xx0(i)-half_win : xx0(i)+half_win,yy0(i)-half_win : yy0(i)+half_win);
    for j = 1:numel(xx1)
        X = fix(xx1(j)) - half_win:ceil(xx1(j)) + half_win;
        Y = fix(yy1(j)) - half_win:ceil(yy1(j)) + half_win;
        im1_win = im1(X,Y);
        [X,Y] = meshgrid(X,Y);
        [Xq,Yq] = meshgrid(xx1(j)-half_win:xx1(j)+half_win,...
                           yy1(j) - half_win:yy1(j) + half_win);
        im1_win_intep = interp2(X,Y,im1_win,Xq,Yq,'linear');
%         im1_win = im1(xx1(j)-half_win : xx1(j)+half_win,yy1(j)-half_win : yy1(j)+half_win);
        xc(i,j) = corr2(im0_win_intep,im1_win_intep);
    end
end
%% One-to-One Matching
while max(xc(:)) > simi_th
    [r,c] = find(xc == max(xc(:)));
    i = r(1);
    j = c(1);
    mx = xc(i,j);
    xc(i,j) = -1;
    next_mx = max(max(xc(i,:),[],2),max(xc(:,j)));
    
    if (1-mx) < (1-next_mx)*dist_th
        crd(i,j) = 1;
    end
    xc(i,:) = -1;
    xc(:,j) = -1;
end
[i,j] = find(crd);
X0 = x0(i);
Y0 = y0(i);
X1 = x1(j);
Y1 = y1(j);

w = abs(X0-X1);
d = abs(Y0-Y1);
n = find(w>20 | d >500);
X0(n) = [];
Y0(n) = [];
X1(n) = [];
Y1(n) = [];
end