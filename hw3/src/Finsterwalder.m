function [p1,p2,p3]=Finsterwalder(a,b,c,j1,j2,j3)
% q: 2D inhomogeneous normalized points(camero coordinate) 2*1
% p: 3D inhomogeneous points(camero coordinate) 3*1
m = zeros(2,1);
n = zeros(2,1);

p1 = [];
p2 = [];
p3 = [];

cos_alpha = dot(j2,j3);
cos_beta = dot(j1,j3);
cos_gamma = dot(j1,j2);

G = c^2*(c^2*(1-cos_beta^2) - b^2*(1-cos_gamma^2));
H = b^2*(b^2 - a^2)*(1-cos_gamma^2) + c^2*(c^2 + 2*a^2)*(1-cos_beta^2)...
    + 2*b^2*c^2*(-1 + cos_alpha*cos_beta*cos_gamma);
I = b^2*(b^2 - c^2)*(1-cos_alpha^2) + a^2*(a^2 + 2*c^2)*(1-cos_beta^2)...
    + 2*a^2*b^2*(-1 + cos_alpha*cos_beta*cos_gamma);
J = a^2*(a^2*(1-cos_beta^2) - b^2*(1-cos_alpha^2));

lambda0 = roots([G H I J]);
lambda_real = [];
for i = 1:numel(lambda0)
    if isreal(lambda0(i))
        lambda_real = [lambda_real,lambda0(i)];
    end
end
%%
for j = 1:numel(lambda_real)
    A = 1 + lambda_real(j);
    B = -cos_alpha;
    C = (b^2 - a^2)/b^2 - lambda_real(j)*c^2/b^2;
    D = -lambda_real(j)*cos_gamma;
    E = (a^2/b^2 + lambda_real(j)*c^2/b^2)*cos_beta;
    F = -a^2/b^2 + lambda_real(j)*(b^2-c^2)/b^2;

    p = sqrt(B^2 - A*C);
    q = sign(B*E - C*D)*sqrt(E^2 - C*F);

    m(1) = (-B + p)/C;
    n(1) = -(E - q)/C;
    m(2) = (-B - p)/C;
    n(2) = -(E + q)/C;

    A = b^2 -m.^2*c^2;
    B = c^2*(cos_beta - n) .* m - b^2*cos_gamma;
    C = -c^2*n.^2 + 2*c^2*n*cos_beta + b^2 - c^2;

    u_large = -sign(B) ./ A .* (abs(B) + sqrt(B.^2 - A .* C));
    u_small = C ./ (A .* u_large);

    u = [];
    for i = 1:2
        if isreal(u_large(i))
            u = [u,[u_large(i);u_small(i)]];
        end
    end
%     u
    if numel(u) ~= 0
        v = u .* m + n;
        u = reshape(u,[],1);
        v = reshape(v,[],1);

        s1 = sqrt(c^2 ./ (1 + u.^2 - 2*u*cos_gamma));
        s2 = u .* s1;
        s3 = v .* s1;

        p1 = s1' .* j1;
        p2 = s2' .* j2;
        p3 = s3' .* j3;
        break
    end
end