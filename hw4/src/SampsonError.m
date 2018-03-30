function delta = SampsonError(x1_inhomo,x2_inhomo,H)
n = size(x1_inhomo,2);
delta = zeros(4,n);
for i = 1:n
    epsilon = [-(x1_inhomo(1,i)*H(2,1)+x1_inhomo(2,i)*H(2,2)+H(2,3))...
               + x2_inhomo(2,i)*(x1_inhomo(1,i)*H(3,1)+x1_inhomo(2,i)*H(3,2)+H(3,3));
                 x1_inhomo(1,i)*H(1,1)+x1_inhomo(2,i)*H(1,2)+H(1,3)....
               - x2_inhomo(1,i)*(x1_inhomo(1,i)*H(3,1)+x1_inhomo(2,i)*H(3,2)+H(3,3))];
    J = [-H(2,1)+x2_inhomo(2,i)*H(3,1), -H(2,2)+x2_inhomo(2,i)*H(3,2),...
         0, x1_inhomo(1,i)*H(3,1)+x1_inhomo(2,i)*H(3,2)+H(3,3);
         H(1,1)-x2_inhomo(1,i)*H(3,1),H(1,2)-x2_inhomo(1,i)*H(3,2),...
         -(x1_inhomo(1,i)*H(3,1)+x1_inhomo(2,i)*H(3,2)+H(3,3)),                                                  0];
    delta(:,i) = -J'*(J*J'\epsilon);
end
end