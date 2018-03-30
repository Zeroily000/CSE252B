function error = SampsonError(x1,x2,F)
n = size(x1,2);
error = zeros(1,n);
for i = 1:n    
    error(i) = (x2(:,i)'*F*x1(:,i))^2  / ...
               ((x2(:,i)'*F(:,1))^2 + (x2(:,i)'*F(:,2))^2 + ...
                (F(1,:)*x1(:,i))^2 + (F(2,:)*x1(:,i))^2);
end
end