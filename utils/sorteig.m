function [V,D] = sorteig(A,order)
if nargin < 2
    order = 'descend';
end

[V1,D1]     = eig(A);
[~,idxsort] = sort(diag(D1),order);
D           = D1(idxsort,idxsort);
V           = V1(:,idxsort);
end