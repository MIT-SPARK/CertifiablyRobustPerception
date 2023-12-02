% Applies the affine map
% SS: R^k -> R^{m x n}
% to a vector u in R^k

function U = applyAffineMap(S,u)

k = length(u);
n = size(S,2);
m = size(S,1)/(k+1);

v = [u(:);1];
U = zeros(m,n);
for i=1:n
    U(:,i) = reshape(S(:,i),[m,k+1])*v;
end
