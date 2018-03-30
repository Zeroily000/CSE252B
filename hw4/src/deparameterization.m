function v_bar = deparameterization(v)
v_norm = sqrt(sum(v.^2));
a = cos(v_norm/2);
b = sin(v_norm/2) ./ v_norm .* v;
v_bar = [a;b];
end