function [z,u,f_est] = recover_solution(X,C,m,k)
X       = X{1};
C       = C{1};
idx     = blkIndices(k+1,m);
X0      = X(idx,idx);
[V,D]   = eig(X0);
[d,idx] = sort(diag(D),'descend');
V       = V(:,idx);
z       = V(:,1);
z       = z/norm(z);

u       = zeros(k,1);
for i = 1:k
    u(i) = trace(X(blkIndices(k+1,m), blkIndices(i,m)));
end

x       = kron([u;1],z);
f_est   = x'*C*x;
end