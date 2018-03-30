function X_PI = TwoViewTriangulation(x1,x2,wu,wv,s,f)
n = size(x1,2);
sigma = Xdeparameterization(s);
D = diag([sigma;0]);
F = expm(skew(wu))*D*expm(skew(wv))';
%% 2D
x1_optimal = zeros(size(x1));
x2_optimal = zeros(size(x2));
for i = 1:n
    %% Fs1
    T1 = [x1(3,i),       0, -x1(1,i)
                0, x1(3,i), -x1(2,i)
                0,       0,  x1(3,i)];
    T2 = [x2(3,i),       0, -x2(1,i)
                0, x2(3,i), -x2(2,i)
                0,       0,  x2(3,i)];
    Fs = (T2' \ F) / T1;
    %% Fs2
    [U,~,V] = svd(Fs);
    e1 = V(:,end);
    e2 = U(:,end);    
    e1 = e1/sqrt(e1(1)^2 + e1(2)^2);
    e2 = e2/sqrt(e2(1)^2 + e2(2)^2);
    R1 = [ e1(1), e1(2), 0
          -e1(2), e1(1), 0
               0,     0, 1];
    R2 = [ e2(1), e2(2), 0
          -e2(2), e2(1), 0
               0,     0, 1];
    Fs = R2*Fs*R1';
    %% min t
    a = Fs(2,2);b=Fs(2,3);c=Fs(3,2);d=Fs(3,3);f1=e1(3);f2=e2(3);
    t = real(roots(f(a,b,c,d,f1,f2)));
    ss = t.^2 ./ (1+f1^2*t.^2) + (c*t+d).^2 ./ ((a*t+b).^2 + f2^2 *(c*t+d).^2);
    [~,id] = min(ss);
    t_min = t(id(end));
    %% closest point
    l1 = [t_min*f1,1,-t_min]';
    l2 = [-f2*(c*t_min+d),a*t_min+b,c*t_min+d]';
    xs1 = [-l1(1)*l1(3),-l1(2)*l1(3),l1(1)^2+l1(2)^2]';
    xs2 = [-l2(1)*l2(3),-l2(2)*l2(3),l2(1)^2+l2(2)^2]';
    %% original coordinate
    x1_optimal(:,i) = (T1\R1')*xs1;
    x2_optimal(:,i) = (T2\R2')*xs2;
end
%% 3D
l2 = F*x1_optimal;
l2_orthogonal = [-l2(2,:) .* x2_optimal(3,:)
                  l2(1,:) .* x2_optimal(3,:)
                  l2(2,:) .* x2_optimal(1,:) - l2(1,:) .* x2_optimal(2,:)];
PI = F2P(wu,wv,s)'*l2_orthogonal;
X_PI = [ PI(4,:).*x1_optimal
        -PI(1,:).*x1_optimal(1,:)-PI(2,:).*x1_optimal(2,:)-PI(3,:).*x1_optimal(3,:)];
end