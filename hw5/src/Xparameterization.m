function v = Xparameterization(v_bar)
v_bar = v_bar ./ sqrt(sum(v_bar.^2));
a = v_bar(1,:);
b = v_bar(2:end,:);
v = 2*acos(a) ./ sin(acos(a)) .* b;
v_norm = sqrt(sum(v.^2));
v = (1 - 2*pi ./ v_norm .* ceil((v_norm-pi)/2/pi)) .* v;
end