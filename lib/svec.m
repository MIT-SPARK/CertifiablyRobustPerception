function v = svec(M)

n = size(M,1);
I = triu(true(n,n),0);
I2 = triu(true(n,n),1);
M(I2') = M(I2')*sqrt(2);
v = M(I');
