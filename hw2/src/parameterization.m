function v = parameterization(v_bar)
a = v_bar(1);
b = v_bar(2:end);
v = 2*acos(a)/sin(acos(a))*b;
if norm(v,2)>pi
    v = (norm(v,2)-2*pi)*v/norm(v,2);
end
end