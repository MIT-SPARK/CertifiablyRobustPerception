% Applies the affine map
% SS: R^k -> R^{m x n}
% to a vector u in R^k

function U = applyAffineMapCell(Scell,u)

k = length(u);

assert(length(Scell)==k+1,'Dimension mismatch.')

U = Scell{end};

for i = 1:k
    U = U + u(i)*Scell{i};
end

U  = full(U);
end
