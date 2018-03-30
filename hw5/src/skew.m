function a_x = skew(a)
a_x = [    0, -a(3),  a(2)
       a(3),     0, -a(1)
      -a(2),  a(1),     0];
end