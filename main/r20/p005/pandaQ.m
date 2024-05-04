y = 0;
n = 1e3;
for i = 2 : n
   y = (2*y + 3)^(1/n);
end
log2(y)