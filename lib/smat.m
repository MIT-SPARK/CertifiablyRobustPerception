function M = smat(v)
dim_v = length(v);
n = (-1 + sqrt(1+8*dim_v))/2;

[si,sj] = meshgrid(1:n,1:n);
si = si(:);
sj = sj(:);

mask = (si >= sj);

si = si(mask);
sj = sj(mask);

vij = (2*n+2 - sj) .* (sj - 1) / 2 + (si - sj) + 1;

mask_diag = (si==sj);
mask_nondiag = ~mask_diag;

si_diag = si(mask_diag);
si_nondiag = si(mask_nondiag);
sj_diag = sj(mask_diag);
sj_nondiag = sj(mask_nondiag);
vij_diag = vij(mask_diag);
vij_nondiag = vij(mask_nondiag);

v_diag = v(vij_diag);
v_nondiag = v(vij_nondiag);
v_nondiag = v_nondiag/sqrt(2);

M = zeros(n,n);
M(sub2ind(size(M),si_nondiag,sj_nondiag)) = v_nondiag;
M = M + M';
M(sub2ind(size(M),si_diag,sj_diag)) = v_diag;



